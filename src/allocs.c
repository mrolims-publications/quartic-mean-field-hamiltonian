#include "../include/allocs.h"
#include <stdio.h>
#include <stdlib.h>

void *xmalloc(size_t n, size_t elem_size) {
    void *ptr = malloc(n * elem_size);
    if (!ptr) {
        fprintf(stderr, "malloc has failed\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void *xcalloc(size_t n, size_t elem_size) {
    void *ptr = calloc(n, elem_size);
    if (!ptr) {
        fprintf(stderr, "calloc has failed\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}
