#ifndef PARSE_ARGS_H
#define PARSE_ARGS_H

#include <stddef.h>

double get_double(size_t argc, char *argv[], size_t pos);
float get_float(size_t argc, char *argv[], size_t pos);
int get_int(size_t argc, char *argv[], size_t pos);
size_t get_size_t(size_t argc, char *argv[], size_t pos);
unsigned int get_uint(size_t argc, char *argv[], size_t pos);

#endif
