#include <stdlib.h>

void memalloc( int *n, void **p )
{
  int msg = posix_memalign(p, 32, (size_t)((size_t)*n *8));
}