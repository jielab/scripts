/* GU CUDA adapter for PhyML. GPL-3.0, derived from PhyML lk.c. */
#include "lk.h"
#include "gu_cuda.h"

int GU_Update_Partial(t_tree *tree, t_edge *b, t_node *d)
{
  t_node *n1=NULL, *n2=NULL;
  phydbl *out=NULL,*a=NULL,*c=NULL,*p=NULL,*q=NULL,*tp=NULL,*tq=NULL;
  int *scale=NULL,*s1=NULL,*s2=NULL,*loc=NULL;
  static int *states=NULL;
  static size_t capacity=0;
  const int n=tree->data->n_pattern, k=tree->mod->ras->n_catg;
  if(tree->mod->ns != 4 || tree->mod->use_m4mod || tree->scaling_method != SCALE_FAST)
    { fprintf(stderr,"[GU PHYML GPU] unsupported model/scaling; aborting\n"); exit(2); }
  Set_All_Partial_Lk(&n1,&n2,&out,&scale,&loc,&p,&tp,&a,&s1,&q,&tq,&c,&s2,d,b,tree);
  if(capacity < (size_t)n) {
    free(states); states=malloc(2*(size_t)n*sizeof(int)); capacity=n;
    if(!states) { fprintf(stderr,"[GU PHYML GPU] host allocation failed\n"); exit(2); }
  }
  for(int i=0;i<n;i++) {
    states[i]=(n1->tax && !n1->c_seq->is_ambigu[i]) ? n1->c_seq->d_state[i] : -1;
    states[n+i]=(n2->tax && !n2->c_seq->is_ambigu[i]) ? n2->c_seq->d_state[i] : -1;
  }
  gu_cuda_partial(n,k,n1->tax,n2->tax,tree->apply_lk_scaling,
                  a,c,p,q,s1,s2,states,states+n,tree->data->wght,out,scale);
  return 1;
}
