#include <stdlib.h>

void memalloc( int *n, int *step, void **p ) 
{
  int msg = posix_memalign(p, 8*(*step), (size_t)((size_t)*n *8));
}