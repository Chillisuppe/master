/***************************************************
* wrapqsort.c
* this is a hack to link gfortran programs to qsort
***************************************************/

#include <stdlib.h>

void qsort_(void* data, int* c, int* s, __compar_fn_t fn)
{
  qsort(data, (size_t)*c, (size_t)*s, fn);
}
