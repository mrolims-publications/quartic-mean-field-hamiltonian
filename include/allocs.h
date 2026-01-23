#ifndef ALLOCS_H
#define ALLOCS_H

#include <stddef.h>

void *xmalloc(size_t n, size_t elem_size);
void *xcalloc(size_t n, size_t elem_size);

#endif
