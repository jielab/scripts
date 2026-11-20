// GU CUDA nucleotide pruning kernel. GPL-3.0; equations/edge cases from PhyML lk.c.
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <cmath>
#include <vector>
#include "gu_cuda.h"
static unsigned long long calls=0, patterns=0;
static void check(cudaError_t e) {
  if(e!=cudaSuccess) { fprintf(stderr,"[GU PHYML GPU] %s\n",cudaGetErrorString(e)); exit(2); }
}
static void report() { fprintf(stderr,"[GU PHYML GPU] partial_kernel_calls=%llu pattern_updates=%llu\n",calls,patterns); }
extern "C" void gu_cuda_check() {
  static bool ready=false; if(ready) return;
  const char *s=getenv("PHYML_GPU_DEVICE"); int device=s?atoi(s):0;
  check(cudaSetDevice(device)); check(cudaFree(0));
  cudaDeviceProp prop; check(cudaGetDeviceProperties(&prop,device));
  fprintf(stderr,"[GU PHYML GPU] backend=cuda-pruning-v1 device=%s precision=double\n",prop.name);
  atexit(report); ready=true;
}
struct Inputs {
  const double *a,*b,*p,*q,*w;
  const int *s1,*s2,*st1,*st2;
};
__global__ void prune(int n,int k,int t1,int t2,int scaling,Inputs in,double *out,int *scale) {
  int site=blockIdx.x*blockDim.x+threadIdx.x;
  if(site>=n) return;
  double *z=out+(size_t)site*k*4;
  if(in.w[site]<=DBL_MIN) { for(int j=0;j<k*4;j++) z[j]=0; scale[site]=0; return; }
  int x=in.st1[site],y=in.st2[site];
  double largest=-DBL_MAX;
  for(int cat=0;cat<k;cat++) {
    const double *a=in.a+((size_t)site*(t1?1:k)+(t1?0:cat))*4;
    const double *b=in.b+((size_t)site*(t2?1:k)+(t2?0:cat))*4;
    bool all_missing=(x<0 && y<0);
    if(all_missing) for(int j=0;j<4;j++) if(a[j]!=1. || b[j]!=1.) all_missing=false;
    for(int i=0;i<4;i++) {
      double value=1.;
      if(!all_missing) {
        const double *p=in.p+cat*16+i*4,*q=in.q+cat*16+i*4;
        double u=0.,v=0.;
        if(x>=0) u=p[x]; else for(int j=0;j<4;j++) u+=p[j]*a[j];
        if(y>=0) v=q[y]; else for(int j=0;j<4;j++) v+=q[j]*b[j];
        value=u*v;
      }
      z[cat*4+i]=value; if(value>largest) largest=value;
    }
  }
  scale[site]=in.s1[site]+in.s2[site];
  if(scaling && largest<0x1p-256) {
    for(int j=0;j<k*4;j++) z[j]*=0x1p256;
    scale[site]+=256;
  }
}
extern "C" void gu_cuda_partial(int n,int k,int t1,int t2,int scaling,
  const double *a,const double *b,const double *p,const double *q,
  const int *s1,const int *s2,const int *st1,const int *st2,
  const double *w,double *out,int *scale) {
  gu_cuda_check();
  if(n<=0 || k<=0) { fprintf(stderr,"Invalid GPU pruning dimensions\n"); exit(2); }
  static std::vector<unsigned char> host;
  static unsigned char *dev=nullptr; static size_t capacity=0;
  size_t na=(size_t)n*(t1?1:k)*4,nb=(size_t)n*(t2?1:k)*4;
  size_t input=(na+nb+32*k+n)*sizeof(double)+4*(size_t)n*sizeof(int);
  size_t output=(size_t)n*k*4*sizeof(double)+(size_t)n*sizeof(int);
  // Ensure the output double buffer remains 8-byte aligned for odd n.
  input=(input+7)&~size_t(7);
  if(input+output>capacity) {
    if(dev) check(cudaFree(dev));
    capacity=input+output; check(cudaMalloc(&dev,capacity)); host.resize(capacity);
  }
  size_t off=0; Inputs in;
  auto pack=[&](const void *src,size_t bytes) {
    size_t pos=off; if(src) memcpy(host.data()+off,src,bytes); else memset(host.data()+off,0,bytes);
    off+=bytes; return dev+pos;
  };
  in.a=(double*)pack(a,na*8); in.b=(double*)pack(b,nb*8);
  in.p=(double*)pack(p,k*16*8); in.q=(double*)pack(q,k*16*8); in.w=(double*)pack(w,n*8);
  in.s1=(int*)pack(s1,n*4); in.s2=(int*)pack(s2,n*4);
  in.st1=(int*)pack(st1,n*4); in.st2=(int*)pack(st2,n*4);
  check(cudaMemcpy(dev,host.data(),input,cudaMemcpyHostToDevice));
  double *result=(double*)(dev+input); int *sc=(int*)(result+(size_t)n*k*4);
  prune<<<(n+127)/128,128>>>(n,k,t1,t2,scaling,in,result,sc); check(cudaGetLastError());
  check(cudaMemcpy(host.data()+input,dev+input,output,cudaMemcpyDeviceToHost));
  memcpy(out,host.data()+input,(size_t)n*k*4*8);
  memcpy(scale,host.data()+input+(size_t)n*k*4*8,(size_t)n*4);
  calls++; patterns+=n;
}
