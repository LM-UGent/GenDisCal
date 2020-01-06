#ifndef DEF_VECOPS
#define DEF_VECOPS

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

typedef double* vec;
typedef int32_t* veci32;
typedef int64_t* veci64;

double* vec_cpy(double* source, size_t len);

typedef double(*vecinfo_f)(double*, size_t);
typedef void(*vecmodify_f)(double*, size_t);
typedef void(*vecipoglobal_f)(double*, double, size_t);
typedef void(*vecipovec_f)(double*, double*, size_t);
typedef double(*vecrelation_f)(double*, double*, size_t);

/* conversions */
int32_t* vec_to_veci32(double* input, size_t len);
int64_t* vec_to_veci64(double* input, size_t len);
double* veci32_to_vec(int32_t* input, size_t len);
double* veci64_to_vec(int64_t* input, size_t len);
int32_t* veci64_to_veci32(int64_t* input, size_t len);
int64_t* veci32_to_veci64(int32_t* input, size_t len);

/* double */
/* extract information */
size_t vec_maxid(double* source, size_t len);
size_t vec_minid(double* source, size_t len);
size_t vec_quantileid(double* source, size_t len, size_t targetid);
size_t vec_medid(double* source, size_t len);
double vec_max(double* source, size_t len);
double vec_min(double* source, size_t len);
double vec_med(double* source, size_t len);
double vec_avg(double* source, size_t len);
double vec_norm(double* source, size_t len);
double vec_sum(double* source, size_t len);

/* in place operations */
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

/* int32_t */
/* extract information */
size_t veci32_minid(int32_t* source, size_t len);
size_t veci32_maxid(int32_t* source, size_t len);
size_t veci32_quantileid(int32_t* source, size_t len, size_t targetid);
size_t veci32_medid(int32_t* source, size_t len);
int32_t veci32_max(int32_t* source, size_t len);
int32_t veci32_min(int32_t* source, size_t len);
int32_t veci32_med(int32_t* source, size_t len);
int32_t veci32_avg(int32_t* source, size_t len);
int32_t veci32_norm(int32_t* source, size_t len);
int32_t veci32_sum(int32_t* source, size_t len);
int64_t veci32_sum64(int32_t* source, size_t len);

/* in place operations */
void veci32_zero(int32_t* target, size_t len);
void veci32_ones(int32_t* target, size_t len);
void veci32_setall(int32_t* target, int32_t val, size_t len);
void veci32_add(int32_t* target, int32_t* B, size_t len);
void veci32_subtract(int32_t* target, int32_t* B, size_t len);
void veci32_add_all(int32_t* target, int32_t B, size_t len);
void veci32_subtract_all(int32_t* target, int32_t B, size_t len);
void veci32_scale(int32_t* target, int32_t factor, size_t len);
void veci32_opposite(int32_t* target, size_t len);
void veci32_dot(int32_t* target, int32_t* B, size_t len);
void veci32_abs(int32_t* target, size_t len);

/* int64_t */
/* extract information */
size_t veci64_minid(int64_t* source, size_t len);
size_t veci64_maxid(int64_t* source, size_t len);
size_t veci64_quantileid(int64_t* source, size_t len, size_t targetid);
size_t veci64_medid(int64_t* source, size_t len);
int64_t veci64_max(int64_t* source, size_t len);
int64_t veci64_min(int64_t* source, size_t len);
int64_t veci64_med(int64_t* source, size_t len);
int64_t veci64_avg(int64_t* source, size_t len);
int64_t veci64_norm(int64_t* source, size_t len);
int64_t veci64_sum(int64_t* source, size_t len);

/* in place operations */
void veci64_zero(int64_t* target, size_t len);
void veci64_ones(int64_t* target, size_t len);
void veci64_setall(int64_t* target, int64_t val, size_t len);
void veci64_add(int64_t* target, int64_t* B, size_t len);
void veci64_subtract(int64_t* target, int64_t* B, size_t len);
void veci64_add_all(int64_t* target, int64_t B, size_t len);
void veci64_subtract_all(int64_t* target, int64_t B, size_t len);
void veci64_scale(int64_t* target, int64_t factor, size_t len);
void veci64_opposite(int64_t* target, size_t len);
void veci64_dot(int64_t* target, int64_t* B, size_t len);
void veci64_abs(int64_t* target, size_t len);


void vec_sort(double* target, size_t len);
void vec_sort_by(void* target, size_t elt_size, double* values, size_t numel);
void vec_sorti64(int64_t* target, size_t len);
void vec_sort_byi64(void* target, size_t elt_size, int64_t* values, size_t numel);
void vec_reorder_byi64(void* target, size_t elt_size, int64_t* new_order, size_t numel);

uint64_t* vec_toranks(double* vec, size_t len);
double vec_variance(double* X, size_t len);
double vec_variance_weighted(double* X, double* W, size_t len);
double vec_covariance(double* X, double* Y, size_t len);
double vec_pearsoncorr(double* X, double* Y, size_t len);
double vec_spearmancorr(double* X, double* Y, size_t len);
double vec_manhattandist(double* A, double* B, size_t len);
double vec_euclidiandist(double* A, double* B, size_t len);


#endif
