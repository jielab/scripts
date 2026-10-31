// SPDX-License-Identifier: GPL-3.0-or-later
// Experimental CUDA port of the IBDmix streaming segment stack.
// Algorithm: PrincetonUniversity/IBDmix, commit 6a9d4bc67dcaf602e65cac0fa6191e1e7a1291ef.
// Frequencies, masks and three per-site LOD values stay on the CPU; sample
// segment stacks and statistics run on CUDA. No change to cohort frequencies.
#include <cuda_runtime.h>
#include <CLI/CLI.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "IBDmix/Mask_Reader.h"
#include "IBDmix/Sample_Mapper.h"
#include "IBDmix/lod_calculator.h"

using U64 = unsigned long long;
static void check(cudaError_t e) {
  if (e != cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));
}
template<class T> struct Device {
  T* p = nullptr;
  explicit Device(size_t n) { check(cudaMalloc(&p, n * sizeof(T))); }
  ~Device() { cudaFree(p); }
  void put(const T* h, size_t n) { check(cudaMemcpy(p,h,n*sizeof(T),cudaMemcpyHostToDevice)); }
  void get(T* h, size_t n) { check(cudaMemcpy(h,p,n*sizeof(T),cudaMemcpyDeviceToHost)); }
};
struct Node { double cumulative, lod; U64 pos; int next; unsigned char bits; };
// Append fixed-size pages instead of reallocating a complete rectangular stack.
// Existing nodes keep their indexes and addresses throughout every extension.
struct NodePages {
  int rows, samples;
  std::vector<Node*> pages;
  Device<Node*> table;
  NodePages(int samples, int initial, int maximum)
    : rows(std::min(initial,2048)), samples(samples), table((maximum+rows-1)/rows) {
    pages.reserve((maximum+rows-1)/rows);
  }
  ~NodePages() { for(Node* page:pages) cudaFree(page); }
  size_t bytes_to(int capacity) const {
    size_t count=(capacity+rows-1)/rows;
    return (count-pages.size())*(size_t)rows*samples*sizeof(Node);
  }
  void grow(int capacity) {
    size_t count=(capacity+rows-1)/rows;
    while(pages.size()<count) {
      Node* page=nullptr; check(cudaMalloc(&page,(size_t)rows*samples*sizeof(Node)));
      pages.push_back(page);
    }
    table.put(pages.data(),pages.size());
  }
};
struct State { int top, start, end, free, used; U64 stats[9]; };
struct Event { U64 start, end, stats[9]; double lod; int row, sample, order; };
struct Site { U64 pos; double lod[3]; unsigned char bits, archaic, selected; };

struct Stack {
  Node** pages; int page_rows; State& s; int sample, samples; int* error;
  __device__ Node& node(int i) { return pages[i/page_rows][(size_t)(i%page_rows)*samples+sample]; }
  __device__ void release(int i) { node(i).next=s.free; s.free=i; s.used--; }
  __device__ int create(U64 pos, double lod, unsigned char bits) {
    int i=s.free;
    if(i<0) { atomicMax(error,1); return -1; }
    s.free=node(i).next;
    s.used++;
    node(i)={lod,lod,pos,-1,bits}; return i;
  }
  __device__ void record(int i) {
    const Node& n=node(i); auto b=n.bits;
    s.stats[0]++; s.stats[1]+=(n.lod>0); s.stats[2]+=(n.lod<0);
    s.stats[3]+=((b&1) && (b&6)); s.stats[4]+=((b&1) && !(b&6));
    s.stats[5]+=(!(b&1) && (b&2)); s.stats[6]+=(!(b&1) && (b&4));
    s.stats[7]+=bool(b&8); s.stats[8]+=bool(b&16);
  }
  __device__ bool add(U64 pos, double lod, unsigned char bits, int row,
                      double threshold, bool exclusive, Event* events,
                      int* count, int capacity) {
    if(s.top<0 && lod<0) return true;
    int todo=create(pos,lod,bits), order=0;
    if(todo<0) return false;
    // The CPU implementation recursively replays the suffix after a maximum.
    // Use a linked work list to preserve exactly that order without recursion.
    while(todo>=0) {
      int i=todo; todo=node(i).next;
      if(s.top<0 && node(i).lod<0) { release(i); continue; }
      node(i).next=s.top; s.top=i;
      if(node(i).next<0) {
        s.start=s.end=i; node(i).cumulative=node(i).lod;
        for(int k=0;k<9;k++) s.stats[k]=0;
        record(i); continue;
      }
      node(i).cumulative=node(node(i).next).cumulative+node(i).lod;
      if(node(i).cumulative>=node(s.end).cumulative) {
        for(int j=i;j!=s.end;j=node(j).next) record(j);
        s.end=i;
        int j=node(i).next;
        while(j!=s.start) { int next=node(j).next; release(j); j=next; }
        node(i).next=s.start;
      }
      if(node(i).cumulative<0) {
        if(node(s.end).cumulative>=threshold) {
          U64 end=node(s.end).pos;
          if(exclusive && s.top!=s.end) {
            int j=s.top; while(node(j).next!=s.end) j=node(j).next;
            if(node(j).lod!=-INFINITY) end=node(j).pos;
          }
          int ix=atomicAdd(count,1);
          if(ix>=capacity) { atomicMax(error,2); return false; }
          Event& e=events[ix]; e.start=node(s.start).pos; e.end=end;
          e.lod=node(s.end).cumulative; e.row=row; e.sample=sample; e.order=order++;
          for(int k=0;k<9;k++) e.stats[k]=s.stats[k];
        }
        // Reverse suffix (top .. end-exclusive), prepend it to pending work.
        int oldtop=s.top, reversed=-1, j=s.top;
        while(j!=s.end) { int next=node(j).next; node(j).next=reversed; reversed=j; j=next; }
        if(reversed>=0) { node(oldtop).next=todo; todo=reversed; }
        j=s.end;
        while(j>=0) { int next=node(j).next; release(j); j=next; }
        s.top=s.start=s.end=-1;
      }
    }
    return true;
  }
};
__global__ void init(Node** pages, int page_rows, State* states, int samples, int capacity) {
  int s=blockIdx.x*blockDim.x+threadIdx.x; if(s>=samples) return;
  states[s]={}; states[s].top=states[s].start=states[s].end=-1; states[s].free=0;
  Stack stack{pages,page_rows,states[s],s,samples,nullptr};
  for(int i=0;i<capacity;i++) stack.node(i).next=(i+1<capacity?i+1:-1);
}
__global__ void extend(Node** pages, int page_rows, State* states, int samples, int old_capacity, int capacity) {
  int s=blockIdx.x*blockDim.x+threadIdx.x; if(s>=samples) return;
  Stack stack{pages,page_rows,states[s],s,samples,nullptr};
  for(int i=old_capacity;i<capacity;i++) stack.node(i).next=(i+1<capacity?i+1:states[s].free);
  states[s].free=old_capacity;
}
__global__ void scan(Node** pages, int page_rows, State* states, const Site* sites,
                     const unsigned char* gt, int samples, int rows,
                     double threshold, bool exclusive, bool purge,
                     Event* events, int* count, int capacity, int* error, int* max_used) {
  int s=blockIdx.x*blockDim.x+threadIdx.x; if(s>=samples) return;
  Stack stack{pages,page_rows,states[s],s,samples,error};
  for(int r=0;r<rows;r++) {
    Site site=sites[r]; unsigned char g=gt[(size_t)r*samples+s], bits=site.bits;
    if(!site.selected && site.archaic==0 && g==2) bits|=16;
    if(!site.selected && site.archaic==2 && g==0) bits|=8;
    double lod=(g==9?0:site.lod[g]);
    if(!stack.add(site.pos,lod,bits,r,threshold,exclusive,events,count,capacity)) return;
  }
  if(purge) stack.add(0,-INFINITY,0,rows,threshold,exclusive,events,count,capacity);
  atomicMax(max_used,states[s].used);
}

static int run(int argc, char** argv) {
  CLI::App app{"Experimental IBDmix CUDA segment caller (CPU frequency/LOD, GPU segment stacks)"};
  std::string genotype_file, output_file="-", sample_file, archaic, mask_file;
  double threshold=3, ae=.01, me=.0025, mp=2;
  int mac=1, batch=2048, nodes=2048, max_nodes=131072, event_capacity=1048576, device=0;
  bool stats=false, inclusive=false, probe=false;
  app.add_option("-g,--genotype",genotype_file);
  app.add_option("-o,--output",output_file);
  app.add_option("-s,--sample",sample_file);
  app.add_option("-n,--archaic",archaic);
  app.add_option("-r,--mask",mask_file);
  app.add_option("-d,--LOD-threshold",threshold)->check(CLI::NonNegativeNumber);
  app.add_option("-m,--minor-allele-count-threshold",mac)->check(CLI::NonNegativeNumber);
  app.add_option("-a,--archaic-error",ae)->check(CLI::Range(0.0,0.499999));
  app.add_option("-e,--modern-error-max",me)->check(CLI::Range(0.0,0.499999));
  app.add_option("-c,--modern-error-proportion",mp)->check(CLI::NonNegativeNumber);
  app.add_flag("-t,--more-stats",stats);
  app.add_flag("-i,--inclusive-end",inclusive);
  app.add_flag("--gpu-check",probe,"Check CUDA device without opening analysis data");
  app.add_option("--gpu-device",device)->check(CLI::NonNegativeNumber);
  app.add_option("--gpu-batch-sites",batch)->check(CLI::Range(1,65536));
  app.add_option("--gpu-nodes",nodes,"Initial stack nodes per selected individual; grows before a batch if needed")->check(CLI::Range(4,1048576));
  app.add_option("--gpu-max-nodes",max_nodes,"Maximum nodes per individual; exceeding the bound fails explicitly")->check(CLI::Range(4,1048576));
  app.add_option("--gpu-events",event_capacity)->check(CLI::Range(1,16777216));
  CLI11_PARSE(app,argc,argv);
  check(cudaSetDevice(device)); cudaDeviceProp prop; check(cudaGetDeviceProperties(&prop,device));
  std::cerr<<"[IBDMIX GPU] backend=cuda-segment-v1 device="<<prop.name<<" double_precision=1\n";
  if(probe) return 0;
  if(genotype_file.empty()) throw std::runtime_error("--genotype is required");
  if(nodes>max_nodes) throw std::runtime_error("--gpu-nodes exceeds --gpu-max-nodes");
  if(!std::isfinite(threshold+ae+me+mp)) throw std::runtime_error("Parameters must be finite");
  auto begin=std::chrono::steady_clock::now();
  std::ifstream input(genotype_file), sample, mask;
  if(!input) throw std::runtime_error("Cannot open genotype file");
  if(!sample_file.empty()) { sample.open(sample_file); if(!sample) throw std::runtime_error("Cannot open sample file"); }
  if(!mask_file.empty()) { mask.open(mask_file); if(!mask) throw std::runtime_error("Cannot open mask file"); }
  std::string line; if(!std::getline(input,line)) throw std::runtime_error("Empty genotype input");
  std::istringstream header(line); Sample_Mapper mapping;
  int n=mapping.initialize(header,sample,archaic);
  if(n<1) throw std::runtime_error("No modern samples selected");
  int columns=0; { std::istringstream h(line); std::string word; while(h>>word) columns++; columns-=4; }
  const auto& names=mapping.getSamples();
  std::vector<int> indexes(n);
  for(int i=0;i<n;i++) indexes[i]=mapping.getSample(i);
  if(std::find(indexes.begin(),indexes.end(),mapping.getArchaicIndex())!=indexes.end())
    throw std::runtime_error("Archaic sample must not be included in modern sample list");
  { auto copy=names; std::sort(copy.begin(),copy.end()); if(std::adjacent_find(copy.begin(),copy.end())!=copy.end()) throw std::runtime_error("Duplicate modern sample names"); }
  Mask_Reader masks(mask_file.empty()?nullptr:&mask); LodCalculator calculator(ae,me,mp,1e-200);
  size_t free_mem,total_mem; check(cudaMemGetInfo(&free_mem,&total_mem));
  NodePages dn(n,nodes,max_nodes);
  size_t required=dn.bytes_to(nodes)+(size_t)n*batch+batch*sizeof(Site)+(size_t)n*sizeof(State)+(size_t)event_capacity*sizeof(Event)+3*sizeof(int);
  if(required>free_mem*0.8) throw std::runtime_error("GPU working memory exceeds 80% of free VRAM; reduce --gpu-nodes/--gpu-batch-sites or use CPU. Do not split the frequency cohort to work around memory limits.");
  std::cerr<<"[IBDMIX GPU] samples="<<n<<" batch_sites="<<batch<<" nodes_per_sample="<<nodes<<" working_MiB="<<required/1048576<<"\n";
  dn.grow(nodes); Device<State> ds(n); Device<Site> dsite(batch);
  Device<unsigned char> dg((size_t)n*batch); Device<Event> de(event_capacity); Device<int> dc(1),derr(1),du(1);
  init<<<(n+127)/128,128>>>(dn.table.p,dn.rows,ds.p,n,nodes); check(cudaGetLastError()); check(cudaDeviceSynchronize());
  std::ofstream out;
  if(output_file!="-") { out.open(output_file); if(!out) throw std::runtime_error("Cannot open output file"); }
  std::ostream& output=(output_file=="-"?std::cout:out);
  output<<"ID\tchrom\tstart\tend\tslod";
  if(stats) output<<"\tsites\tpositive_lods\tnegative_lods\tmask_and_maf\tin_mask\tmaf_low\tmaf_high\trec_2_0\trec_0_2";
  output<<'\n';
  std::vector<Site> sites(batch); std::vector<unsigned char> genotypes((size_t)batch*n);
  std::vector<Event> events; std::string chromosome; U64 previous=0, total_sites=0, total_events=0;
  double gpu_seconds=0; auto last_log=begin; int max_used=0;
  for(;;) {
    int rows=0;
    while(rows<batch && std::getline(input,line)) {
      if(line.empty()) throw std::runtime_error("Blank genotype record");
      size_t tabs[4], at=0;
      for(int j=0;j<4;j++) { tabs[j]=line.find('\t',at); if(tabs[j]==std::string::npos) throw std::runtime_error("Malformed genotype record"); at=tabs[j]+1; }
      std::string chrom=line.substr(0,tabs[0]);
      if(chromosome.empty()) chromosome=chrom;
      if(chrom!=chromosome) throw std::runtime_error("GPU caller accepts one chromosome per call");
      size_t used=0; std::string position=line.substr(tabs[0]+1,tabs[1]-tabs[0]-1);
      U64 pos=std::stoull(position,&used);
      if(used!=position.size() || !pos || pos<previous) throw std::runtime_error("Invalid or unsorted genotype position"); previous=pos;
      const char* data=line.data()+at; size_t width=line.size()-at;
      if(width && data[width-1]=='\r') width--;
      if(width && data[width-1]=='\t') width--; // upstream generate_gt trailing tab
      if(width!=(size_t)columns*2-1) throw std::runtime_error("Genotype columns differ from header (expected canonical tab-separated 0/1/2/9)");
      for(int j=0;j<columns;j++) {
        char c=data[2*j]; if(c!='0' && c!='1' && c!='2' && c!='9') throw std::runtime_error("Invalid genotype code");
        if(j+1<columns && data[2*j+1]!='\t') throw std::runtime_error("Invalid genotype separator");
      }
      U64 alt=0, total=0;
      for(int j=0;j<n;j++) { unsigned char g=data[2*indexes[j]]-'0'; genotypes[(size_t)rows*n+j]=g; if(g!=9) { alt+=g; total+=2; } }
      Site site{}; site.pos=pos; site.archaic=data[2*mapping.getArchaicIndex()]-'0';
      if(masks.in_mask(chromosome,pos)) site.bits|=1;
      if(alt<=(U64)mac) site.bits|=2;
      if(total-alt<=(U64)mac) site.bits|=4;
      site.selected=(site.bits==0);
      calculator.update_lod_cache(site.archaic+'0',total?(double)alt/total:0,site.selected);
      for(int j=0;j<3;j++) site.lod[j]=calculator.calculate_lod(j+'0');
      sites[rows++]=site;
    }
    if(input.bad()) throw std::runtime_error("Genotype read failed");
    bool final=input.eof(); int zero=0;
    // At most one new node per input site. Replay uses existing nodes. Reserve
    // this upper bound before the kernel, preserving every in-progress stack.
    int needed=max_used+rows+(final?1:0);
    if(needed>nodes) {
      int capacity=std::min(((needed+dn.rows-1)/dn.rows)*dn.rows,max_nodes);
      if(capacity<needed) throw std::runtime_error("GPU stack capacity exceeded; increase --gpu-max-nodes (no valid completed output)");
      check(cudaMemGetInfo(&free_mem,&total_mem));
      size_t additional=dn.bytes_to(capacity);
      if(additional>free_mem*0.8) throw std::runtime_error("GPU stack growth exceeds free VRAM: additional_MiB="+std::to_string(additional/1048576)+" free_MiB="+std::to_string(free_mem/1048576)+"; use CPU without changing the frequency cohort");
      dn.grow(capacity);
      extend<<<(n+127)/128,128>>>(dn.table.p,dn.rows,ds.p,n,nodes,capacity);
      check(cudaGetLastError()); check(cudaDeviceSynchronize());
      std::cerr<<"[IBDMIX GPU] grow_nodes="<<nodes<<"->"<<capacity<<" added_MiB="<<additional/1048576<<"\n"; nodes=capacity;
    }
    dc.put(&zero,1); derr.put(&zero,1); du.put(&zero,1); dsite.put(sites.data(),rows); dg.put(genotypes.data(),(size_t)rows*n);
    auto gpu_start=std::chrono::steady_clock::now();
    scan<<<(n+127)/128,128>>>(dn.table.p,dn.rows,ds.p,dsite.p,dg.p,n,rows,threshold,!inclusive,final,de.p,dc.p,event_capacity,derr.p,du.p);
    check(cudaGetLastError()); check(cudaDeviceSynchronize());
    gpu_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-gpu_start).count();
    int error=0,count=0; derr.get(&error,1); dc.get(&count,1); du.get(&max_used,1);
    if(error) throw std::runtime_error(error==1?"GPU stack capacity exceeded; rerun with larger --gpu-nodes (no valid completed output)":"GPU event capacity exceeded; increase --gpu-events or reduce --gpu-batch-sites (no valid completed output)");
    events.resize(count); de.get(events.data(),count);
    std::sort(events.begin(),events.end(),[](const Event& a,const Event& b) { if(a.row!=b.row) return a.row<b.row; if(a.sample!=b.sample) return a.sample<b.sample; return a.order<b.order; });
    for(const Event& e:events) {
      output<<names[e.sample]<<'\t'<<chromosome<<'\t'<<e.start<<'\t'<<e.end<<'\t'<<e.lod;
      if(stats) for(auto v:e.stats) output<<'\t'<<v;
      output<<'\n';
    }
    total_sites+=rows; total_events+=count; auto now=std::chrono::steady_clock::now();
    if(final || std::chrono::duration<double>(now-last_log).count()>=30) {
      std::cerr<<"[IBDMIX GPU] sites="<<total_sites<<" segments="<<total_events<<" elapsed_s="<<std::chrono::duration<double>(now-begin).count()<<" gpu_kernel_s="<<gpu_seconds<<"\n"; last_log=now;
    }
    if(final) break;
  }
  output.flush(); if(!output) throw std::runtime_error("Output write failed");
  return 0;
}
int main(int argc,char** argv) {
  try { return run(argc,argv); }
  catch(const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; return 2; }
}
