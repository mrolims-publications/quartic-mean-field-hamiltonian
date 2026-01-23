#include "../include/parse_args.h"
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

double get_double(size_t argc, char *argv[], size_t pos) {
    if (pos >= argc) {
        fprintf(stderr, "Missing argument %zu\n", pos);
        exit(EXIT_FAILURE);
    }
    char *endptr;
    errno = 0;
    double val = strtod(argv[pos], &endptr);
    if (errno == ERANGE || *endptr != '\0') {
        fprintf(stderr, "Invalid double at argv[%zu]: '%s'\n", pos, argv[pos]);
        exit(EXIT_FAILURE);
    }
    return val;
}

float get_float(size_t argc, char *argv[], size_t pos) {
    if (pos >= argc) {
        fprintf(stderr, "Missing argument %zu\n", pos);
        exit(EXIT_FAILURE);
    }
    char *endptr;
    errno = 0;
    float val = strtof(argv[pos], &endptr);
    if (errno == ERANGE || *endptr != '\0') {
        fprintf(stderr, "Invalid float at argv[%zu]: '%s'\n", pos, argv[pos]);
        exit(EXIT_FAILURE);
    }
    return val;
}

int get_int(size_t argc, char *argv[], size_t pos) {
    if (pos >= argc) {
        fprintf(stderr, "Missing argument %zu\n", pos);
        exit(EXIT_FAILURE);
    }
    char *endptr;
    errno = 0;
    long val = strtol(argv[pos], &endptr, 10);
    if (errno == ERANGE || *endptr != '\0' || val > INT_MAX || val < INT_MIN) {
        fprintf(stderr, "Invalid int at argv[%zu]: '%s'\n", pos, argv[pos]);
        exit(EXIT_FAILURE);
    }
    return (int)val;
}

unsigned int get_uint(size_t argc, char *argv[], size_t pos) {
    if (pos >= argc) {
        fprintf(stderr, "Missing argument %zu\n", pos);
        exit(EXIT_FAILURE);
    }
    char *endptr;
    errno = 0;
    long val = strtol(argv[pos], &endptr, 10);
    if (errno == ERANGE || *endptr != '\0' || val > INT_MAX || val < INT_MIN) {
        fprintf(stderr, "Invalid int at argv[%zu]: '%s'\n", pos, argv[pos]);
        exit(EXIT_FAILURE);
    }
    return (unsigned int)val;
}

size_t get_size_t(size_t argc, char *argv[], size_t pos) {
    if (pos >= argc) {
        fprintf(stderr, "Missing argument %zu\n", pos);
        exit(EXIT_FAILURE);
    }
    char *endptr;
    errno = 0;
    unsigned long long val = strtoull(argv[pos], &endptr, 10);
    if (errno == ERANGE || *endptr != '\0') {
        fprintf(stderr, "Invalid size_t at argv[%zu]: '%s'\n", pos, argv[pos]);
        exit(EXIT_FAILURE);
    }
    return (size_t)val;
}
