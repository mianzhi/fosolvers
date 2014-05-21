//---------------------------------------------------------------------------- best with 100 columns

#include <stdlib.h>

/**
 * node
 */
typedef struct node node;
struct node{
  node* p; /// parent
  node* c[8]; /// children
  int v; /// payload
};

node* init();
void purge(node **pself);
void clear(node *self);
void split(node *self);
void splitId(node *root, long i);
void trim(node *self);
void trimId(node *root, long i);
void grow(node *self, int lvl);
int nBinDigits(long i);
void getPos(long i, double p[3]);
int neibId(long i, int d);
int nLf(node *self);
void getNode(node *self, long i, node **ptar);
void getNxtLf(node *root, long *i, node **pnxt);
void markLf(node *root);
void genConn(node *root, long id[], int lvl[], int neib[][6]);

