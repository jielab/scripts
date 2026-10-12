// CUDA kernels for TRACE gamma emissions and pairwise transition posteriors.
// Equations follow tracehmm/trace.py and trace_utils.pyx; MIT (see LICENSE).
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
static std::string error;
static unsigned long long emission_calls=0,transition_calls=0;
static void check(cudaError_t e) { if(e!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(e)); }
static void report() { fprintf(stderr,"[GU TRACE GPU] emission_kernel_calls=%llu transition_kernel_calls=%llu\n",emission_calls,transition_calls); }
extern "C" const char *gu_trace_error() { return error.c_str(); }
extern "C" int gu_trace_check(int device) {
  try {
    check(cudaSetDevice(device)); check(cudaFree(0));
    cudaDeviceProp p; check(cudaGetDeviceProperties(&p,device));
    fprintf(stderr,"[GU TRACE GPU] backend=cuda-emission-transition-v1 device=%s precision=double\n",p.name);
    atexit(report); return 0;
  } catch(const std::exception &e) { error=e.what(); return 1; }
}
static double *buf=nullptr; static size_t capacity=0;
static void allocate(size_t n) {
  if(n>capacity) { if(buf) check(cudaFree(buf)); buf=nullptr; capacity=0; check(cudaMalloc(&buf,n*sizeof(double))); capacity=n; }
}
__global__ void emissions(int n,const double *x,const double *mask,double *out,
  double a1,double b1,double ab1,double lg1,double a2,double b2,double ab2,double lg2,
  double cp,double left,double right,double logs) {
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n) return;
  double value=0.;
  if(mask[i]!=0.) {
    if(cp>=0. && x[i]>=cp) value=(a1-1.)*log(x[i])-b1*x[i]+ab1-lg1+right;
    else value=(a2-1.)*log(x[i])-b2*x[i]+ab2-lg2+left;
  }
  out[i]=value+logs;
}
extern "C" int gu_trace_emissions(int n,const double *x,const double *mask,double *out,
  double a1,double b1,double ab1,double lg1,double a2,double b2,double ab2,double lg2,
  double cp,double left,double right,double logs) {
  try {
    allocate(3*(size_t)n);
    check(cudaMemcpy(buf,x,n*8,cudaMemcpyHostToDevice));
    check(cudaMemcpy(buf+n,mask,n*8,cudaMemcpyHostToDevice));
    emissions<<<(n+255)/256,256>>>(n,buf,buf+n,buf+2*n,a1,b1,ab1,lg1,a2,b2,ab2,lg2,cp,left,right,logs);
    check(cudaGetLastError()); check(cudaMemcpy(out,buf+2*n,n*8,cudaMemcpyDeviceToHost));
    emission_calls++; return 0;
  } catch(const std::exception &e) { error=e.what(); return 1; }
}
__global__ void transitions(int n,const double *a,const double *b,const double *e,
  double lp,double lq,double l1p,double l1q,double *o1,double *o2) {
  int i=blockIdx.x*blockDim.x+threadIdx.x; if(i>=n-1) return; int j=i+1;
  double v[4]={a[i]+l1p+b[j]+e[j],a[i]+lp+b[n+j]+e[n+j],
               a[n+i]+lq+b[j]+e[j],a[n+i]+l1q+b[n+j]+e[n+j]};
  double m=-1e32; for(int k=0;k<4;k++) if(v[k]>m) m=v[k];
  double c=0.; for(int k=0;k<4;k++) c+=exp(v[k]-m);
  double norm=m+log(c); o1[i]=v[1]-norm; o2[i]=v[2]-norm;
}
extern "C" int gu_trace_transitions(int n,const double *a,const double *b,const double *e,
  double lp,double lq,double l1p,double l1q,double *o1,double *o2) {
  try {
    allocate(8*(size_t)n);
    check(cudaMemcpy(buf,a,2*(size_t)n*8,cudaMemcpyHostToDevice));
    check(cudaMemcpy(buf+2*n,b,2*(size_t)n*8,cudaMemcpyHostToDevice));
    check(cudaMemcpy(buf+4*n,e,2*(size_t)n*8,cudaMemcpyHostToDevice));
    if(n>1) {
      transitions<<<(n+254)/256,256>>>(n,buf,buf+2*n,buf+4*n,lp,lq,l1p,l1q,buf+6*n,buf+7*n);
      check(cudaGetLastError());
      check(cudaMemcpy(o1,buf+6*n,(n-1)*8,cudaMemcpyDeviceToHost));
      check(cudaMemcpy(o2,buf+7*n,(n-1)*8,cudaMemcpyDeviceToHost));
    }
    transition_calls++; return 0;
  } catch(const std::exception &e) { error=e.what(); return 1; }
}
