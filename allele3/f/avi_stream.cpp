// Bounded-memory merged-Atlas VCF scanner. One chromosome per worker.
// Build: c++ -O3 -std=c++17 avi_stream.cpp -lz -o avi_stream
#include <zlib.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
struct Reader {
 gzFile f; char buf[65536];
 explicit Reader(const string&p){f=gzopen(p.c_str(),"rb");if(!f)throw runtime_error("Cannot open "+p);gzbuffer(f,4*1024*1024);}
 bool line(string&s){s.clear();while(gzgets(f,buf,sizeof(buf))){s+=buf;if(!s.empty()&&s.back()=='\n'){s.pop_back();if(!s.empty()&&s.back()=='\r')s.pop_back();return true;}}int e=0;const char*m=gzerror(f,&e);if(e!=Z_OK&&e!=Z_STREAM_END)throw runtime_error(m);return !s.empty();}
 ~Reader(){gzclose(f);}
};
vector<string> split(const string&s,char sep){vector<string>r;size_t a=0,b;while((b=s.find(sep,a))!=string::npos){r.push_back(s.substr(a,b-a));a=b+1;}r.push_back(s.substr(a));return r;}
string tag(const string&s,const string&key){string find=key+"=";size_t a=s.find(find);while(a!=string::npos&&a>0&&s[a-1]!=';')a=s.find(find,a+1);if(a==string::npos)return ".";a+=find.size();size_t b=s.find(';',a);return s.substr(a,b==string::npos?b:b-a);}
vector<double> scores(const string&s,size_t n){auto p=split(s,',');vector<double>v(n,NAN);if(s==".")return v;if(p.size()!=n)throw runtime_error("Number=A score cardinality differs from ALT");for(size_t i=0;i<n;i++){if(p[i]=="."||p[i].empty())continue;size_t end;v[i]=stod(p[i],&end);if(end!=p[i].size()||!isfinite(v[i]))throw runtime_error("Invalid score");}return v;}
string canonical(string c){if(c.rfind("chr",0)==0)c=c.substr(3);if(c=="MT")c="M";return "chr"+c;}
struct Acc {
 long size,start=-1,end=0,npos=0,ncomplete=0,nsub=0,nhigh=0,nlow=0,nadj=0;
 double sum=0,sum2=0,sm=0,sm2=0,srange=0,ssd=0,grad=0,grad2=0;
 Acc(long w):size(w){}
 void add(long pos,const vector<double>&v,bool adjacent,double difference){
  if(start<0){start=((pos-1)/size)*size;end=start+size;}
  npos++;if(v.size()==3)ncomplete++;nsub+=v.size();
  double s=0,s2=0;for(double x:v){s+=x;s2+=x*x;}double mean=s/v.size();
  sum+=s;sum2+=s2;sm+=mean;sm2+=mean*mean;
  srange+=*max_element(v.begin(),v.end())-*min_element(v.begin(),v.end());ssd+=sqrt(max(0.,s2/v.size()-mean*mean));
  if(adjacent){nadj++;grad+=abs(difference);grad2+=difference*difference;}
 }
};
void val(ostream&o,double x){if(isfinite(x))o<<x;else o<<"NA";}
void emit(ostream&o,const string&chrom,const Acc&a,long length){if(!a.npos)return;long end=length>0?min(a.end,length):a.end;double ns=a.nsub,np=a.npos;
 o<<chrom<<'\t'<<a.start<<'\t'<<end<<'\t'<<a.size<<'\t'<<a.npos<<'\t'<<a.ncomplete<<'\t'<<a.nsub<<'\t'<<np/(end-a.start)<<'\t'<<a.sum/ns<<'\t'<<sqrt(max(0.,a.sum2/ns-pow(a.sum/ns,2)))<<'\t'<<a.sm/np<<'\t'<<sqrt(max(0.,a.sm2/np-pow(a.sm/np,2)))<<'\t'<<a.srange/np<<'\t'<<a.ssd/np<<'\t'<<a.nadj<<'\t';val(o,a.nadj?a.grad/a.nadj:NAN);o<<'\t';val(o,a.nadj?sqrt(a.grad2/a.nadj):NAN);o<<'\t'<<a.nhigh/ns<<'\t'<<a.nlow/ns<<'\n';}
int main(int argc,char**argv){try{
 map<string,string>a;for(int i=1;i<argc;i+=2){if(i+1>=argc)throw runtime_error("All options need a value");a[argv[i]]=argv[i+1];}
 for(auto k:{"--input","--out","--chrom"})if(!a.count(k))throw runtime_error(string("Required ")+k);
 string chrom=canonical(a["--chrom"]),field=a.count("--field")?a["--field"]:"AG_AVI_PHRED";
 bool query=a.count("--targets");vector<long>targets;size_t target_i=0;
 if(query){Reader r(a["--targets"]);string s;while(r.line(s)){if(s.empty()||s[0]=='#')continue;try{size_t n;long p=stol(s,&n);if(p<1||n!=s.size())throw runtime_error("Invalid target");targets.push_back(p);}catch(const exception&){throw runtime_error("Targets must contain one positive 1-based position per line");}}sort(targets.begin(),targets.end());targets.erase(unique(targets.begin(),targets.end()),targets.end());}
 vector<Acc>windows;if(!query)for(auto s:split(a.count("--windows")?a["--windows"]:"1000,10000,100000,1000000",',')){long w=stol(s);if(w<1)throw runtime_error("Positive window sizes required");windows.emplace_back(w);}
 double hi=a.count("--high")?stod(a["--high"]):20.,lo=a.count("--low")?stod(a["--low"]):5.;
 ofstream out(a["--out"]);if(!out)throw runtime_error("Cannot write output");out<<setprecision(12);
 if(query)out<<"chrom\tpos\tref\talts\tavi_raw\tavi_phred\tsplice\n";
 else out<<"chrom\tstart\tend\twindow_bp\tn_positions\tn_complete_positions\tn_substitutions\tcoverage_fraction\tmean_avi\tsd_avi\tmean_site_avi\tsd_site_mean\tmean_alt_range\tmean_alt_sd\tn_adjacent_pairs\tmean_abs_gradient\trms_gradient\thigh_fraction\tlow_fraction\n";
 Reader in(a["--input"]);string s;long lastpos=-1,previous_score_pos=-1,length=0;double prev=NAN;long parsed=0,matched=0,unscored=0;bool header=false;
 while(in.line(s)){
  if(s.empty())continue;
  if(s[0]=='#'){
   if(s.rfind("##reference=",0)==0&&(s.find("GRCh37")!=string::npos||s.find("hg19")!=string::npos))throw runtime_error("Expected GRCh38 Atlas");
   if(s.rfind("##contig=",0)==0){size_t q=s.find("ID=");if(q!=string::npos){q+=3;size_t e=s.find_first_of(",>",q);if(canonical(s.substr(q,e-q))==chrom){size_t p=s.find("length=");if(p!=string::npos){length=stol(s.substr(p+7));if(length<1)throw runtime_error("Invalid contig length");}}}}
   if(s.rfind("#CHROM",0)==0)header=true;continue;
  }
  if(!header)throw runtime_error("VCF #CHROM header missing");
  size_t x=0;vector<string>v;for(int k=0;k<7;k++){size_t y=s.find('\t',x);if(y==string::npos)throw runtime_error("VCF needs 8 columns");v.push_back(s.substr(x,y-x));x=y+1;}
  if(canonical(v[0])!=chrom)continue;
  size_t consumed;long pos=stol(v[1],&consumed);if(consumed!=v[1].size()||pos<1)throw runtime_error("Invalid 1-based position");if(length>0&&pos>length)throw runtime_error("Position exceeds declared contig length");if(pos<=lastpos)throw runtime_error("Merged Atlas must contain one sorted record per position");lastpos=pos;parsed++;
  if(query){while(target_i<targets.size()&&targets[target_i]<pos)target_i++;if(target_i==targets.size())break;if(targets[target_i]!=pos)continue;}
  auto alts=split(v[4],',');bool snv=v[3].size()==1&&string("ACGT").find(v[3])!=string::npos;
  for(auto alt:alts)snv=snv&&alt.size()==1&&string("ACGT").find(alt)!=string::npos&&alt!=v[3];
  if(!snv)continue;
  auto distinct=alts;sort(distinct.begin(),distinct.end());if(unique(distinct.begin(),distinct.end())!=distinct.end())throw runtime_error("Duplicate ALT bases in merged Atlas");
  string info=s.substr(x,s.find('\t',x)==string::npos?string::npos:s.find('\t',x)-x);
  if(query){auto raw=tag(info,"AG_AVI_RAW"),phred=tag(info,"AG_AVI_PHRED");scores(raw,alts.size());scores(phred,alts.size());auto splice=tag(info,"AG_SPLICE");scores(splice,alts.size());out<<chrom<<'\t'<<pos<<'\t'<<v[3]<<'\t'<<v[4]<<'\t'<<raw<<'\t'<<phred<<'\t'<<splice<<'\n';matched++;continue;}
  auto values=scores(tag(info,field),alts.size());vector<double>good;for(double t:values)if(isfinite(t))good.push_back(t);
  if(good.empty()){unscored++;continue;}
  double mean=0;for(double t:good)mean+=t;mean/=good.size();
  for(auto&w:windows){if(w.start>=0&&pos-1>=w.end){emit(out,chrom,w,length);w=Acc(w.size);}bool adjacent=pos==previous_score_pos+1&&isfinite(prev)&&((pos-2)/w.size==(pos-1)/w.size);w.add(pos,good,adjacent,mean-prev);for(double t:good){w.nhigh+=t>=hi;w.nlow+=t<=lo;}}
  previous_score_pos=pos;prev=mean;matched++;
 }
 if(!query)for(auto&w:windows)emit(out,chrom,w,length);
 if(!out.good())throw runtime_error("Output write failed");
 cerr<<"AVI "<<chrom<<" parsed="<<parsed<<" matched="<<matched<<" unscored="<<unscored<<'\n';return 0;
 }catch(const exception&e){cerr<<"avi_stream: "<<e.what()<<'\n';return 2;}}
