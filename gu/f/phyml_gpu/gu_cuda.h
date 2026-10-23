#ifndef GU_PHYML_CUDA_H
#define GU_PHYML_CUDA_H
#ifdef __cplusplus
extern "C" {
#endif
void gu_cuda_check(void);
void gu_cuda_partial(int n,int k,int tip1,int tip2,int scaling,
                     const double *a,const double *b,const double *p,const double *q,
                     const int *s1,const int *s2,const int *st1,const int *st2,
                     const double *w,double *out,int *scale);
#ifdef __cplusplus
}
#endif
#endif
