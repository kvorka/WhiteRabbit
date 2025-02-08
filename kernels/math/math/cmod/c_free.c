#include <stdlib.h>

void memfree( void **p )
{
  free(*p);
}