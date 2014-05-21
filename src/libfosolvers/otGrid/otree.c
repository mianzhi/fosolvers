//---------------------------------------------------------------------------- best with 100 columns

#include <stdlib.h>
#include <stdio.h>
#include "otree.h"

int const BPD=3; /// bit per digit
int const I_HEAD=04; /// index prefix
int const N_FACE=6; /// number of faces

/**
 * initialize node
 */
node* init()
{
  node *ptr=NULL;
  ptr=malloc(sizeof(node));
  ptr->p=NULL;
  for(int i=0;i<8;i++){
    ptr->c[i]=NULL;
  }
  ptr->v=0;
  return ptr;
}

/**
 * recursively purge node
 */
void purge(node **pself)
{
  if((*pself)->c[0]!=NULL){
    for(int i=0;i<8;i++){
      purge(&((*pself)->c[i]));
    }
  }
  free(*pself);
  *pself=NULL;
  return;
}

/**
 * recursively purge node without reset root pointer
 */
void clear(node *self)
{
  if(self->c[0]!=NULL){
    for(int i=0;i<8;i++){
      purge(&(self->c[i]));
    }
  }
  free(self);
  return;
}

/**
 * create children
 */
void split(node *self)
{
  if(self->c[0]==NULL){
    for(int i=0;i<8;i++){
      self->c[i]=init();
      self->c[i]->p=self;
    }
  }else{
    printf("[E] split(): node has children\n");
  }
  return;
}

/**
 * create children with octal id
 */
void splitId(node *root, long i)
{
  node *ptr=NULL;
  getNode(root,i,&ptr);
  split(ptr);
  return;
}

/**
 * remove children
 */
void trim(node *self)
{
  if(self->c[0]!=NULL){
    for(int i=0;i<8;i++){
      purge(&(self->c[i]));
    }
  }else{
    printf("[W] trim(): try trimming leaf\n");
  }
  return;
}

/**
 * remove children with octal id
 */
void trimId(node *root, long i)
{
  node *ptr=NULL;
  getNode(root,i,&ptr);
  trim(ptr);
  return;
}

/**
 * grow the node by given levels
 */
void grow(node *self, int lvl)
{
  if(lvl>0){
    split(self);
    if(lvl>1){
      for(int i=0;i<8;i++){
        grow(self->c[i],lvl-1);
      }
    }
  }else{
    printf("[W] grow(): invalid lvl=%d\n",lvl);
  }
  return;
}

/**
 * number of binary digit in an integer
 */
int nBinDigits(long i)
{
  int n=0;
  for(;i!=0;i>>=1){
    n++;
  }
  return n;
}


/**
 * find position of a octal id
 */
void getPos(long i, double p[BPD])
{
  int n=nBinDigits(i);
  double h=0.25;
  p[0]=p[1]=p[2]=0.5;
  for(int j=n-2*BPD;j>=0;j-=BPD){
    for(int k=0;k<3;k++){
      int l=j+k;
      p[k]+=h*(double)(2*((i&(1<<l))>>l)-1);
    }
    h*=0.5;
  }
  return;
}

/**
 * octal id of neighbor cell
 */
int neibId(long i, int d)
{
  int n=nBinDigits(i);
  int c=(d&6)>>1; /// x, y, z
  int s=(d&1); /// -, +
  if(d>=N_FACE||d<0){
    printf("[E] neibId(): invalid direction d=%d\n",d);
  }
  int f=0;
  int j;
  for(j=c;j<=n-1-BPD;j+=BPD){
    if((i&(1<<j))>>j!=s){
      i^=1<<j;
      f=1;
      break;
    }
  }
  if(f){
    for(j-=BPD;j>0;j-=BPD){
      if((i&(1<<j))>>j==s){
        i^=1<<j;
      }
    }
  }else{
    i=I_HEAD;
  }
  return i;
}

/**
 * total number of leaves
 */
int nLf(node *self)
{
  int n=0;
  if(self->c[0]==NULL){
    n=1;
  }else{
    for(int i=0;i<8;i++){
      n+=nLf(self->c[i]);
    }
  }
  return n;
}

/**
 * get node pointer by index
 */
void getNode(node *self, long i, node **ptar)
{
  if(i==I_HEAD||self->c[0]==NULL){
    *ptar=self;
  }else{
    int n=nBinDigits(i);
    int k=(i>>(n-BPD*2))%8;
    i=(1<<(n-BPD-1))+i%(1<<(n-BPD*2));
    getNode(self->c[k],i,ptar);
  }
  return;
}

/**
 * get next leaf node pointer
 */
void getNxtLf(node *root, long *pi, node **pnxt)
{
  if(*pi==I_HEAD){
    while((*pnxt)->c[0]!=NULL){
      (*pi)<<=BPD;
      *pnxt=(*pnxt)->c[0];
    }
  }else{
    (*pi)++;
    while(*pi%8==0){
      *pi>>=BPD;
    }
    int n=nBinDigits(*pi);
    if((*pi)>>(n-BPD)==I_HEAD+1){
      *pi=I_HEAD;
      *pnxt=NULL;
    }else{
      getNode(root,*pi,pnxt);
      while((*pnxt)->c[0]!=NULL){
        (*pi)<<=BPD;
        *pnxt=(*pnxt)->c[0];
      }
    }
  }
  return;
}

/**
 * mark leaf in the payload with non-zero index
 */
void markLf(node *root)
{
  node *ptr=root;
  long i=I_HEAD;
  int j=0;
  while(ptr!=NULL){
    if(ptr->c[0]==NULL){
      j++;
      ptr->v=j;
    }
    getNxtLf(root,&i,&ptr);
  }
}

/**
 * generate connectivity table
 */
void genConn(node *root, long id[], int lvl[], int neib[][N_FACE])
{
  node *ptr=root;
  long i=I_HEAD;
  int j=0;
  markLf(root);
  while(ptr!=NULL){
    if(ptr->c[0]==NULL){
      id[j]=i;
      lvl[j]=nBinDigits(i)/BPD-1;
      for(int k=0;k<N_FACE;k++){
        node *pn=NULL;
        getNode(root,neibId(i,k),&pn);
        neib[j][k]=pn->v;
      }
      j++;
    }
    getNxtLf(root,&i,&ptr);
  }
}

