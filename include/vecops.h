#ifndef DEF_VECOPS
#define DEF_VECOPS

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

double* vec_cpy(double* source, size_t len);

// extract information
double vec_max(double* source, size_t len);
double vec_min(double* source, size_t len);
double vec_med(double* source, size_t len);
double vec_avg(double* source, size_t len);
double vec_norm(double* source, size_t len);
double vec_sum(double* source, size_t len);

// in place operations
void vec_zero(double* target, size_t len);
void vec_ones(double* target, size_t len);
void vec_setall(double* target, double val, size_t len);
void vec_add(double* target, double* B, size_t len);
void vec_subtract(double* target, double* B, size_t len);
void vec_add_all(double* target, double B, size_t len);
void vec_subtract_all(double* target, double B, size_t len);
void vec_scale(double* target, double factor, size_t len);
void vec_opposite(double* target, size_t len);
void vec_inverse(double* target, size_t len);
void vec_dot(double* target, double* B, size_t len);
void vec_abs(double* target, size_t len);
void vec_sort(double* target, size_t len);


#endif
