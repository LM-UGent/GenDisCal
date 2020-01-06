/*
MIT License

Copyright (c) 2019 Gleb Goussarov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "vecops.h"
#include <string.h>

double* vec_cpy(double* source, size_t len) {
    double* result;
    size_t i = 0;
    result = malloc(sizeof(double)*len);
    for (i = 0;i < len;i++) result[i] = source[i];
    return result;
}

/* conversions */
int32_t* vec_to_veci32(double* input, size_t len){
    int32_t* result;
    size_t i;
    result = malloc(sizeof(int32_t)*len);
    for (i = 0;i < len;i++)
        result[i] = (int32_t)(input[i]);
    return result;
}
int64_t* vec_to_veci64(double* input, size_t len){
    int64_t* result;
    size_t i;
    result = malloc(sizeof(int64_t)*len);
    for (i = 0;i < len;i++)
        result[i] = (int64_t)(input[i]);
    return result;
}
double* veci32_to_vec(int32_t* input, size_t len){
    double* result;
    size_t i;
    result = malloc(sizeof(double)*len);
    for (i = 0;i < len;i++)
        result[i] = (double)(input[i]);
    return result;
}
double* veci64_to_vec(int64_t* input, size_t len){
    double* result;
    size_t i;
    result = malloc(sizeof(double)*len);
    for (i = 0;i < len;i++)
        result[i] = (double)(input[i]);
    return result;
}
int32_t* veci64_to_veci32(int64_t* input, size_t len){
    int32_t* result;
    size_t i;
    result = malloc(sizeof(int32_t)*len);
    for (i = 0;i < len;i++)
        result[i] = (int32_t)(input[i]);
    return result;
}
int64_t* veci32_to_veci64(int32_t* input, size_t len){
    int64_t* result;
    size_t i;
    result = malloc(sizeof(int64_t)*len);
    for (i = 0;i < len;i++)
        result[i] = (int64_t)(input[i]);
    return result;
}

size_t vec_maxid(double* source, size_t len) {
    size_t i;
    size_t id;
    double value = 0.0;
    id = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++) {
        if (value < source[i]) {
            value = source[i];
            id = i;
        }
    }
    return id;
}
size_t vec_minid(double* source, size_t len) {
    size_t i;
    size_t id;
    double value = 0.0;
    id = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++) {
        if (value > source[i]) {
            value = source[i];
            id = i;
        }
    }
    return id;
}
size_t vec_quantileid(double* source, size_t len, size_t targid) {
    /* uses a variation of quick-sort, where only the relevant portion is sorted */
    size_t* neworder;
    size_t* oldorder;
    size_t* tmpptr;
    size_t result;
    double pivot_val;
    double minval, maxval;
    size_t i, i0, imax, nextleft, nextright;
    int sorted;

    if (targid >= len) return vec_maxid(source, len);
    sorted = 1;
    for (i = 1;i < len;i++) {
        if (source[i] < source[i - 1]) {
            sorted = 0;
            break;
        }
    }
    if (sorted) return targid;

    oldorder = malloc(sizeof(size_t)*len);
    neworder = malloc(sizeof(size_t)*len);

    minval = source[0];
    maxval = source[0];
    oldorder[0] = 0;
    sorted = 0;
    for (i = 1;i < len;i++) {
        if (minval > source[i]) minval = source[i];
        if (maxval < source[i]) maxval = source[i];
        oldorder[i] = i;
    }
    pivot_val = minval + (maxval - minval)*((double)(targid) / (double)(len));

    i0 = 0;
    imax = len;
    /* iterate until either the leftmost or the rightmost element matches the target quantile */
    while (imax > targid + 1 && i0 < targid && !sorted) {
        nextleft = i0;
        nextright = imax - 1;
        /* At the end of the following loop, all elements smaller than pivot_val are on the left of nextleft */
        sorted = 1;
        for (i = i0;i < imax;i++) {
            if (sorted && i > i0 && source[oldorder[i]] > source[oldorder[i - 1]]) sorted = 0;
            if (source[oldorder[i]] < pivot_val) {
                neworder[nextleft] = oldorder[i];
                nextleft++;
            }
            else {
                neworder[nextright] = oldorder[i];
                nextright--;
            }
        }
        if (nextleft < targid) i0 = nextleft+1;
        else if (nextleft == targid) i0 = nextleft;
        else imax = nextleft;
        tmpptr = oldorder;
        oldorder = neworder;
        neworder = tmpptr;
        pivot_val = source[oldorder[i0]];
    }
    if (sorted) result = oldorder[targid];
    else {
        if (i0 == targid) {
            /* leftmost element matches quantile - find the smallest element */
            result = oldorder[i0];
            for (i = i0 + 1;i < imax;i++) {
                if (source[oldorder[i]] < source[result])result = oldorder[i];
            }
        }
        else if (imax == targid + 1) {
            result = oldorder[i0];
            /* rightmost element matches quantile - find the biggest element */
            for (i = i0 + 1;i < imax;i++) {
                if (source[oldorder[i]] > source[result])result = oldorder[i];
            }
        }
        else {
            /* multiple elements seem to match - take the central one */
            result = oldorder[(i0 + imax) / 2];
        }
    }
    free(neworder);
    free(oldorder);
    return result;
}
size_t vec_medid(double* source, size_t len) {
    return vec_quantileid(source, len, len / 2);
}

double vec_max(double * source, size_t len) {
    size_t i;
    double value = 0.0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++)
        value = value > source[i] ? value : source[i];
    return value;
}
double vec_min(double * source, size_t len) {
    size_t i;
    double value = 0.0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++)
        value = value < source[i] ? value : source[i];
    return value;
}
double vec_med(double* source, size_t len) {
    return source[vec_medid(source,len)];
}
double vec_avg(double* source, size_t len) {
    size_t i;
    double value = 0.0;
    for (i = 0;i < len;i++)
        value += source[i];
    value /= (double)len;
    return value;
}
double vec_norm(double* source, size_t len) {
    size_t i;
    double value = 0.0;
    for (i = 0;i < len;i++) value += source[i] * source[i];
    value = sqrt(value);
    return value;
}
double vec_sum(double* source, size_t len) {
    size_t i;
    double value = 0.0;
    for (i = 0;i < len;i++)
        value += source[i];
    return value;
}

void vec_zero(double* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = 0.0;
}
void vec_ones(double* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = 1.0;
}
void vec_setall(double* target, double val, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = val;
}
void vec_add(double* target, double * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] += B[i];
}
void vec_add_all(double* target, double B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] += B;
}
void vec_subtract(double* target, double * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] -= B[i];
}
void vec_subtract_all(double* target, double B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] -= B;
}
void vec_scale(double* target, double factor, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] *= factor;
}
void vec_opposite(double* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = -target[i];
}
void vec_inverse(double* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = 1.0 / target[i];
}
void vec_dot(double* target, double * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] *= B[i];
}
void vec_abs(double* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = target[i] > 0.0 ? target[i] : -target[i];
}

/* same instructions for i32 vectors, with the exception of inverse, which does not make sense in the context of integers */
int32_t* veci32_cpy(int32_t* source, size_t len) {
    int32_t* result;
    size_t i = 0;
    result = malloc(sizeof(int32_t)*len);
    for (i = 0;i < len;i++) result[i] = source[i];
    return result;
}

size_t veci32_maxid(int32_t* source, size_t len) {
    size_t i;
    size_t id;
    int32_t value = 0;
    id = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++) {
        if (value > source[i]) {
            value = source[i];
            id = i;
        }
    }
    return id;
}
size_t veci32_minid(int32_t* source, size_t len) {
    size_t i;
    size_t id;
    int32_t value = 0;
    id = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++) {
        if (value < source[i]) {
            value = source[i];
            id = i;
        }
    }
    return id;
}
size_t veci32_quantileid(int32_t* source, size_t len, size_t targid) {
    /* uses a variation of quick-sort, where only the relevant portion is sorted */
    size_t* neworder;
    size_t* oldorder;
    size_t* tmpptr;
    size_t result;
    int32_t pivot_val;
    int32_t minval, maxval;
    size_t i, i0, imax, nextleft, nextright;
    int sorted;

    if (targid >= len) return veci32_maxid(source, len);
    sorted = 1;
    for (i = 1;i < len;i++) {
        if (source[i] < source[i - 1]) {
            sorted = 0;
            break;
        }
    }
    if (sorted) return targid;

    oldorder = malloc(sizeof(size_t)*len);
    neworder = malloc(sizeof(size_t)*len);

    minval = source[0];
    maxval = source[0];
    oldorder[0] = 0;
    sorted = 0;
    for (i = 1;i < len;i++) {
        if (minval > source[i]) minval = source[i];
        if (maxval < source[i]) maxval = source[i];
        oldorder[i] = i;
    }
    pivot_val = minval + (maxval - minval)*((double)(targid) / (double)(len));

    i0 = 0;
    imax = len;
    /* iterate until either the leftmost or the rightmost element matches the target quantile */
    while (imax > targid + 1 && i0 < targid && !sorted) {
        nextleft = i0;
        nextright = imax - 1;
        /* At the end of the following loop, all elements smaller than pivot_val are on the left of nextleft */
        sorted = 1;
        for (i = i0;i < imax;i++) {
            if (sorted && i > i0 && source[oldorder[i]] > source[oldorder[i - 1]]) sorted = 0;
            if (source[oldorder[i]] < pivot_val) {
                neworder[nextleft] = oldorder[i];
                nextleft++;
            }
            else {
                neworder[nextright] = oldorder[i];
                nextright--;
            }
        }
        if (nextleft < targid) i0 = nextleft + 1;
        else if (nextleft == targid) i0 = nextleft;
        else imax = nextleft;
        tmpptr = oldorder;
        oldorder = neworder;
        neworder = tmpptr;
        pivot_val = source[oldorder[i0]];
    }
    if (sorted) result = oldorder[targid];
    else {
        if (i0 == targid) {
            /* leftmost element matches quantile - find the smallest element */
            result = oldorder[i0];
            for (i = i0 + 1;i < imax;i++) {
                if (source[oldorder[i]] < source[result])result = oldorder[i];
            }
        }
        else if (imax == targid + 1) {
            result = oldorder[i0];
            /* rightmost element matches quantile - find the biggest element */
            for (i = i0 + 1;i < imax;i++) {
                if (source[oldorder[i]] > source[result])result = oldorder[i];
            }
        }
        else {
            /* multiple elements seem to match - take the central one */
            result = oldorder[(i0 + imax) / 2];
        }
    }
    free(neworder);
    free(oldorder);
    return result;
}
size_t veci32_medid(int32_t* source, size_t len) {
    return veci32_quantileid(source, len, len / 2);
}
int32_t veci32_max(int32_t * source, size_t len) {
    size_t i;
    int32_t value = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++)
        value = value > source[i] ? value : source[i];
    return value;
}
int32_t veci32_min(int32_t * source, size_t len) {
    size_t i;
    int32_t value = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++)
        value = value < source[i] ? value : source[i];
    return value;
}
int32_t veci32_med(int32_t* source, size_t len) {
    return source[veci32_medid(source, len)];
}
int32_t veci32_avg(int32_t* source, size_t len) {
    size_t i;
    int32_t value = 0;
    for (i = 0;i < len;i++)
        value += source[i];
    value /= (int32_t)len;
    return value;
}
int32_t veci32_norm(int32_t* source, size_t len) {
    size_t i;
    double value = 0;
    for (i = 0;i < len;i++) value += source[i] * source[i];
    value = sqrt(value);
    return (int32_t)value;
}
int32_t veci32_sum(int32_t* source, size_t len) {
    size_t i;
    int32_t value = 0;
    for (i = 0;i < len;i++)
        value += source[i];
    return value;
}
int64_t veci32_sum64(int32_t* source, size_t len) {
    size_t i;
    int64_t value = 0;
    for (i = 0;i < len;i++)
        value += source[i];
    return value;
}

void veci32_zero(int32_t* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = 0;
}
void veci32_ones(int32_t* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = 1;
}
void veci32_setall(int32_t* target, int32_t val, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = val;
}
void veci32_add(int32_t* target, int32_t * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] += B[i];
}
void veci32_add_all(int32_t* target, int32_t B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] += B;
}
void veci32_subtract(int32_t* target, int32_t * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] -= B[i];
}
void veci32_subtract_all(int32_t* target, int32_t B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] -= B;
}
void veci32_scale(int32_t* target, int32_t factor, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] *= factor;
}
void veci32_opposite(int32_t* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = -target[i];
}
void veci32_dot(int32_t* target, int32_t * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] *= B[i];
}
void veci32_abs(int32_t* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = target[i] > 0 ? target[i] : -target[i];
}

/* same instructions for i64 vectors, with the exception of inverse, which does not make sense in the context of integers */
int64_t* veci64_cpy(int64_t* source, size_t len) {
    int64_t* result;
    size_t i = 0;
    result = malloc(sizeof(int64_t)*len);
    for (i = 0;i < len;i++) result[i] = source[i];
    return result;
}

size_t veci64_maxid(int64_t* source, size_t len) {
    size_t i;
    size_t id;
    int64_t value = 0;
    id = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++) {
        if (value > source[i]) {
            value = source[i];
            id = i;
        }
    }
    return id;
}
size_t veci64_minid(int64_t* source, size_t len) {
    size_t i;
    size_t id;
    int64_t value = 0;
    id = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++) {
        if (value < source[i]) {
            value = source[i];
            id = i;
        }
    }
    return id;
}
size_t veci64_quantileid(int64_t* source, size_t len, size_t targid) {
    /* uses a variation of quick-sort, where only the relevant portion is sorted */
    size_t* neworder;
    size_t* oldorder;
    size_t* tmpptr;
    size_t result;
    int64_t pivot_val;
    int64_t minval, maxval;
    size_t i, i0, imax, nextleft, nextright;
    int sorted;

    if (targid >= len) return veci64_maxid(source, len);
    sorted = 1;
    for (i = 1;i < len;i++) {
        if (source[i] < source[i - 1]) {
            sorted = 0;
            break;
        }
    }
    if (sorted) return targid;

    oldorder = malloc(sizeof(size_t)*len);
    neworder = malloc(sizeof(size_t)*len);

    minval = source[0];
    maxval = source[0];
    oldorder[0] = 0;
    sorted = 0;
    for (i = 1;i < len;i++) {
        if (minval > source[i]) minval = source[i];
        if (maxval < source[i]) maxval = source[i];
        oldorder[i] = i;
    }
    pivot_val = minval + (maxval - minval)*((double)(targid) / (double)(len));

    i0 = 0;
    imax = len;
    /* iterate until either the leftmost or the rightmost element matches the target quantile */
    while (imax > targid + 1 && i0 < targid && !sorted) {
        nextleft = i0;
        nextright = imax - 1;
        /* At the end of the following loop, all elements smaller than pivot_val are on the left of nextleft */
        sorted = 1;
        for (i = i0;i < imax;i++) {
            if (sorted && i > i0 && source[oldorder[i]] > source[oldorder[i - 1]]) sorted = 0;
            if (source[oldorder[i]] < pivot_val) {
                neworder[nextleft] = oldorder[i];
                nextleft++;
            }
            else {
                neworder[nextright] = oldorder[i];
                nextright--;
            }
        }
        if (nextleft < targid) i0 = nextleft + 1;
        else if (nextleft == targid) i0 = nextleft;
        else imax = nextleft;
        tmpptr = oldorder;
        oldorder = neworder;
        neworder = tmpptr;
        pivot_val = source[oldorder[i0]];
    }
    if (sorted) result = oldorder[targid];
    else {
        if (i0 == targid) {
            /* leftmost element matches quantile - find the smallest element */
            result = oldorder[i0];
            for (i = i0 + 1;i < imax;i++) {
                if (source[oldorder[i]] < source[result])result = oldorder[i];
            }
        }
        else if (imax == targid + 1) {
            result = oldorder[i0];
            /* rightmost element matches quantile - find the biggest element */
            for (i = i0 + 1;i < imax;i++) {
                if (source[oldorder[i]] > source[result])result = oldorder[i];
            }
        }
        else {
            /* multiple elements seem to match - take the central one */
            result = oldorder[(i0 + imax) / 2];
        }
    }
    free(neworder);
    free(oldorder);
    return result;
}
size_t veci64_medid(int64_t* source, size_t len) {
    return veci64_quantileid(source, len, len / 2);
}

int64_t veci64_max(int64_t * source, size_t len) {
    size_t i;
    int64_t value = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++)
        value = value > source[i] ? value : source[i];
    return value;
}
int64_t veci64_min(int64_t * source, size_t len) {
    size_t i;
    int64_t value = 0;
    if (len>0) value = source[0];
    for (i = 1;i < len;i++)
        value = value < source[i] ? value : source[i];
    return value;
}
int64_t veci64_med(int64_t* source, size_t len) {
    return source[veci64_medid(source, len)];
}
int64_t veci64_avg(int64_t* source, size_t len) {
    size_t i;
    int64_t value = 0;
    for (i = 0;i < len;i++)
        value += source[i];
    value /= (int64_t)len;
    return value;
}
int64_t veci64_norm(int64_t* source, size_t len) {
    size_t i;
    double value = 0.0;
    for (i = 0;i < len;i++) value += source[i] * source[i];
    value = sqrt(value);
    return (int64_t)value;
}
int64_t veci64_sum(int64_t* source, size_t len) {
    size_t i;
    int64_t value = 0;
    for (i = 0;i < len;i++)
        value += source[i];
    return value;
}

void veci64_zero(int64_t* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = 0;
}
void veci64_ones(int64_t* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = 1;
}
void veci64_setall(int64_t* target, int64_t val, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = val;
}
void veci64_add(int64_t* target, int64_t * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] += B[i];
}
void veci64_add_all(int64_t* target, int64_t B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] += B;
}
void veci64_subtract(int64_t* target, int64_t * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] -= B[i];
}
void veci64_subtract_all(int64_t* target, int64_t B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] -= B;
}
void veci64_scale(int64_t* target, int64_t factor, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] *= factor;
}
void veci64_opposite(int64_t* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = -target[i];
}
void veci64_dot(int64_t* target, int64_t * B, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] *= B[i];
}
void veci64_abs(int64_t* target, size_t len) {
    size_t i = 0;
    for (i = 0;i < len;i++) target[i] = target[i] > 0 ? target[i] : -target[i];
}

#define arrswapids(arr,index1,index2,tmpvar)     tmpvar = arr[index1]; arr[index1] = arr[index2]; arr[index2] = tmpvar
static void _vec_quicksort(double* target, size_t len) {
    size_t pivot;
    double pivotvalue;
    double tmp;
    size_t i;

    if (len < 2) return;

    pivot = 1;
    pivotvalue = target[0];
    for (i = 1;i < len; i++) {
        if (target[i] <= pivotvalue) {
            arrswapids(target, i, pivot, tmp);
            pivot++;
        }
    }
    if (pivot == len) {
        pivot--;
        arrswapids(target, 0, pivot, tmp);
        /* doing this prevents an infinite loop */
    }
    _vec_quicksort(target, pivot);
    _vec_quicksort(target+pivot, len-pivot);
}
static void _vec_heapify(double* target, size_t len, size_t node) {
    size_t newroot = node;
    size_t left, right;
    double tmp;
    left = (node << 1) + 1;
    right = (node << 1) + 2;
    if (left < len && target[left] > target[newroot])
        newroot = left;
    if (right < len && target[right] > target[newroot])
        newroot = right;
    if (newroot != node) {
        arrswapids(target, node, newroot, tmp);
        _vec_heapify(target, len, newroot);
    }
}
static void _vec_sort__wdepth(double* target, size_t len, size_t depth, size_t maxdepth) {
    size_t pivot;
    double pivotvalue;
    double tmp;
    size_t i,j,k;

    if (len < 16) {
        for (i = 1;i < len;i++) {
            j = i;
            while (j > 0 && target[j] < target[j - 1])
                j--;
            tmp = target[i];
            for (k = i;k > j;k--) {
                target[k] = target[k - 1];
            }
            target[j] = tmp;
        }
    }
    else if(depth<maxdepth){
        pivot = 1;
        pivotvalue = target[0];
        for (i = 1;i < len; i++) {
            if (target[i] <= pivotvalue) {
                arrswapids(target, i, pivot, tmp);
                pivot++;
            }
        }
        if (pivot == len) {
            pivot--;
            arrswapids(target, 0, pivot, tmp);
        }
        _vec_quicksort(target, pivot);
        _vec_quicksort(target + pivot, len - pivot);
    }
    else {
        for (i = len / 2 - 1; i >= 0; i--) {
            _vec_heapify(target, len, i);
        }
        for (i = len - 1; i >= 0; i--) {
            arrswapids(target, 0, i, tmp);
            _vec_heapify(target, i, 0);
        }
    }
}
void vec_sort(double* target, size_t len) {
    _vec_sort__wdepth(target, len, 0, 32);
}
static void _vec_quicksort_by(void* target, size_t elt_size, double* values, size_t len) {
    size_t pivot;
    double pivotvalue;
    double tmp;
    char part_c;
    char* target_c;
    int sorted;
    size_t i,j,k,k0,k1;
    if (len < 2) return;

    target_c = target;
    sorted = 1;
    for (i = 1;i < len; i++) {
        if (values[i] < values[i - 1]) {
            sorted = 0;
            break;
        }
    }
    if (sorted) return;

    if (len == 2) {
        if (values[0] > values[1]) {
            arrswapids(values, 0, 1, tmp);
            k0 = 0;
            k1 = elt_size;
            for (j = 0;j < elt_size;j++) {
                arrswapids(target_c, k0 + j, k1 + j, part_c);
            }
        }
        return;
    }

    pivot = (len + 1) / 2 - 1;
    pivotvalue = values[pivot];
    tmp = 0;

    if (values[pivot] < values[0]) {
        arrswapids(values, 0, pivot, tmp);
        k0 = 0; k1 = pivot*elt_size;
        for (j = 0;j < elt_size;j++) { arrswapids(target_c, k0 + j, k1 + j, part_c); }
    }
    if (values[len - 1] < values[0]) {
        arrswapids(values, 0, len - 1, tmp);
        k0 = 0; k1 = (len - 1)*elt_size;
        for (j = 0;j < elt_size;j++) { arrswapids(target_c, k0 + j, k1 + j, part_c); }
    }
    if (values[pivot] < values[len - 1]) {
        arrswapids(values, pivot, len - 1, tmp);
        k0 = pivot*elt_size; k1 = (len - 1)*elt_size;
        for (j = 0;j < elt_size;j++) { arrswapids(target_c, k0 + j, k1 + j, part_c); }
    }
    pivotvalue = values[len - 1];

    i = 0;
    j = len - 1;
    while (values[i] < pivotvalue)i++;
    while (j>0 && values[j] >= pivotvalue)j--;
    if (i < j) {
        while (i < j) {
            arrswapids(values, i, j, tmp);
            k0 = i*elt_size; k1 = j*elt_size;
            for (k = 0;k < elt_size;k++) { arrswapids(target_c, k0 + k, k1 + k, part_c); }
            i++;
            if (j>0)j--;
            while (values[i] < pivotvalue)i++;
            while (j>0 && values[j] > pivotvalue)j--;
        }
    }
    else if (i == 0) {
        i++;
    }
    pivot = i;
    _vec_quicksort_by(target, elt_size, values, pivot);
    _vec_quicksort_by((void*)(((char*)target)+elt_size*pivot), elt_size, values + pivot, len - pivot);
}
void vec_sort_by(void* target, size_t elt_size, double* values, size_t numel) {
    _vec_quicksort_by(target, elt_size, values, numel);
}

void _vec_quicksorti64(int64_t* target, size_t len) {
    size_t pivot;
    int64_t pivotvalue;
    int64_t tmp;
    size_t i,j;

    if (len < 2) return;
    if (len == 2) {
        if (target[0] > target[1]) {
            arrswapids(target, 0, 1, tmp);
        }
        return;
    }

    pivot = (len+1)/2-1;
    pivotvalue = target[pivot];
    tmp = 0;

    if (target[pivot] < target[0]) {
        arrswapids(target, 0, pivot, tmp);
    }
    if (target[len - 1] < target[0]) {
        arrswapids(target, 0, len - 1, tmp);
    }
    if (target[pivot] < target[len - 1]) {
        arrswapids(target, pivot, len - 1, tmp);
    }
    pivotvalue = target[len - 1];

    i = 0;
    j = len - 1;
    while (target[i] < pivotvalue)i++;
    while (j>0 && target[j] >= pivotvalue)j--;
    if (i < j) {
        while (i < j) {
            arrswapids(target, i, j, tmp);
            i++;
            if(j>0)j--;
            while (target[i] < pivotvalue)i++;
            while (j>0 && target[j] > pivotvalue)j--;
        }
    }
    pivot = i;
    _vec_quicksorti64(target, pivot);
    _vec_quicksorti64(target + pivot, len - pivot);
}
size_t _vec_quicksorti64_nofollowup(int64_t* target, size_t len) {
    size_t pivot;
    int64_t pivotvalue;
    int64_t tmp;
    size_t i, j;

    if (len < 2) return 0;
    if (len == 2) {
        if (target[0] > target[1]) {
            arrswapids(target, 0, 1, tmp);
        }
        return 0;
    }

    pivot = (len + 1) / 2 - 1;
    pivotvalue = target[pivot];
    tmp = 0;

    if (target[pivot] < target[0]) {
        arrswapids(target, 0, pivot, tmp);
    }
    if (target[len - 1] < target[0]) {
        arrswapids(target, 0, len - 1, tmp);
    }
    if (target[pivot] < target[len - 1]) {
        arrswapids(target, pivot, len - 1, tmp);
    }
    pivotvalue = target[len - 1];

    i = 0;
    j = len - 1;
    while (target[i] < pivotvalue)i++;
    while (j>0 && target[j] >= pivotvalue)j--;
    if (i < j) {
        while (i < j) {
            arrswapids(target, i, j, tmp);
            i++;
            if (j>0)j--;
            while (target[i] < pivotvalue)i++;
            while (j>0 && target[j] > pivotvalue)j--;
        }
    }
    pivot = i;
    return pivot;
}

void _vec_heapifyi64(int64_t* target, size_t len, size_t node) {
    size_t newroot = node;
    size_t left, right;
    int64_t tmp;
    left = (node << 1) + 1;
    right = (node << 1) + 2;
    if (left < len && target[left] > target[newroot])
        newroot = left;
    if (right < len && target[right] > target[newroot])
        newroot = right;
    if (newroot != node) {
        arrswapids(target, node, newroot, tmp);
        _vec_heapifyi64(target, len, newroot);
    }
}
void _vec_sort__wdepthi64(int64_t* target, size_t len, size_t depth, size_t maxdepth) {
    size_t pivot;
    int64_t tmp;
    size_t i, j, k;

    if (len < 8) {
        /* insertion sort */
        for (i = 1;i < len;i++) {
            j = i;
            while (j > 0 && target[i] < target[j - 1])
                j--;
            tmp = target[i];
            for (k = i;k > j;k--) {
                target[k] = target[k - 1];
            }
            target[j] = tmp;
        }
    }
    else if (depth<maxdepth) {
        pivot = _vec_quicksorti64_nofollowup(target, len);
        _vec_sort__wdepthi64(target, pivot, depth+1, maxdepth);
        _vec_sort__wdepthi64(target + pivot, len - pivot, depth + 1, maxdepth);
    }
    else {
        /* note i is unsigned, so if it ever gets bigger than its starting value,
           that means it's reached 0 */
        for (i = len / 2 - 1; i < len; i--) {
            _vec_heapifyi64(target, len, i);
        }
        for (i = len - 1; i < len; i--) {
            arrswapids(target, 0, i, tmp);
            _vec_heapifyi64(target, i, 0);
        }
    }
}
void vec_sorti64(int64_t* target, size_t len) {
    _vec_sort__wdepthi64(target, len, 0, 32);
}
void vec_reorder_byi64(void* target, size_t elt_size, int64_t* new_order, size_t numel) {
    char* tmp;
    size_t i;
    tmp = (char*)malloc(elt_size*numel);
    memcpy(tmp, target, elt_size*numel);
    for (i = 0;i < numel;i++) {
        memcpy(((char*)target) + elt_size*i, tmp + elt_size*new_order[i], elt_size);
    }
    free(tmp);
}


static inline void _vec_insertsortbyvec(void* target, size_t elt_size, double* values, size_t count) {
    size_t i, j, insertat;
    double curvalue;
    int64_t* neworder;
    neworder = (int64_t*)malloc(sizeof(int64_t)*count);
    for (i = 0;i < count;i++) {
        neworder[i] = (int64_t)i;
    }
    for (i = 1;i < count;i++) {
        insertat = i;
        curvalue = values[i];
        j = i;
        while (j > 0 && values[j - 1] > curvalue) {
            neworder[j] = neworder[j - 1];
            values[j] = values[j - 1];
            j--;
        }
        neworder[j] = i;
        values[j] = curvalue;
    }
    vec_reorder_byi64(target, elt_size, neworder, count);
    free(neworder);
}

static inline void _vec_insertsortbyveci64(void* target, size_t elt_size, int64_t* values, size_t count) {
    size_t i, j, insertat;
    int64_t curvalue;
    int64_t* neworder;
    neworder = (int64_t*)malloc(sizeof(int64_t)*count);
    for (i = 0;i < count;i++) {
        neworder[i] = (int64_t)i;
    }
    for (i = 1;i < count;i++) {
        insertat = i;
        curvalue = values[i];
        j = i;
        while (j > 0 && values[j - 1] > curvalue) {
            neworder[j] = neworder[j - 1];
            values[j] = values[j - 1];
            j--;
        }
        neworder[j] = i;
        values[j] = curvalue;
    }
    vec_reorder_byi64(target, elt_size, neworder, count);
    free(neworder);
}

static inline void _vec_insertsortbyi64(void* target, size_t elt_size, int64_t* values, size_t count) {
    size_t i, j, insertat;
    int64_t curvalue;
    int64_t* neworder;
    neworder = (int64_t*)malloc(sizeof(int64_t)*count);
    for (i = 0;i < count;i++) {
        neworder[i] = (int64_t)i;
    }
    for (i = 1;i < count;i++) {
        insertat = i;
        curvalue = values[i];
        j = i;
        while (j > 0 && values[j] > curvalue) {
            neworder[j] = neworder[j - 1];
            values[j] = values[j - 1];
            j--;
        }
        neworder[j] = i;
        values[j] = curvalue;
    }
    vec_reorder_byi64(target, elt_size, neworder, count);
}
static inline size_t _vec_quicksortbyi64_makepivot(void* target, size_t elt_size, int64_t* values, size_t begin, size_t end) {
    size_t pivot, i, j, k0, k1;
    int64_t pivotvalue, tmp;
    char* target_c;
    char part_c;
    int sorted;
    if ((int64_t)end - (int64_t)begin < 2)return end;
    if ((int64_t)end - (int64_t)begin < 16) {

    }
    pivotvalue = values[begin];
    sorted = 1;
    pivot = begin;
    for (i = begin + 1;i < end; i++) {
        if (values[i] < values[i - 1]) sorted = 0;
        if (values[i] <= pivotvalue) {
            /* Whenever a smaller value is discovered, it takes the place of the pivot. */
            /* This ensures that a different pivot is selected each time this function is called. */
            arrswapids(values, i, pivot, tmp);
            k0 = i*elt_size;
            k1 = pivot*elt_size;
            for (j = 0;j < elt_size;j++) {
                arrswapids(target_c, k0 + j, k1 + j, part_c);
            }
            pivot++;
        }
    }
    if (pivot == begin) pivot++;
    if (pivot == end) pivot--;
    if (sorted) pivot = end;
    /* result: every element with a value <= pivotvalue will be below "pivot" */
    /* pivot==begin+1 <=> the first value was the smallest, but the array is not sorted */
    /* pivot==end-1 <=> the first element was the biggest */
    /* pivot==end <=> the array is sorted */
    return pivot;
}
static void _vec_quicksort_byi64(void* target, size_t elt_size, int64_t* values, size_t len) {
    size_t pivot;
    int64_t pivotvalue;
    int64_t tmp;
    char part_c;
    char* target_c;
    size_t i, j, k0, k1;
    if (len < 2) return;

    target_c = target;
    pivot = 1;
    pivotvalue = values[0];
    for (i = 1;i < len; i++) {
        if (values[i] <= pivotvalue) {
            arrswapids(values, i, pivot, tmp);
            k0 = i*elt_size;
            k1 = pivot*elt_size;
            for (j = 0;j < elt_size;j++) {
                arrswapids(target_c, k0 + j, k1 + j, part_c);
            }
            pivot++;
        }
    }
    if (pivot == len) {
        pivot--;
        arrswapids(values, 0, pivot, tmp);
        k1 = pivot*elt_size;
        for (j = 0;j < elt_size;j++) {
            arrswapids(target_c, j, k1 + j, part_c);
        }
        /* doing this prevents an infinite loop */
    }
    _vec_quicksort_byi64(target, elt_size, values, pivot);
    _vec_quicksort_byi64((void*)(((char*)target) + elt_size*pivot), elt_size, values + pivot, len - pivot);
}
void vec_sort_byi64(void* target, size_t elt_size, int64_t* values, size_t numel) {
    _vec_quicksort_byi64(target, elt_size, values, numel);
}
/*
void vec_quicksort_byi64_parallel(void* target, size_t elt_size, int64_t* values, size_t numel, size_t nthreads) {
    int* p_activethreads;
    size_t* nregions;
    size_t* nalloc_rb;
    size_t** regionsbuffer;
    size_t t;
    *p_activethreads = 1;
    nregions = (size_t*) calloc(nthreads, sizeof(size_t));
    regionsbuffer = (size_t**) malloc(sizeof(size_t*)*nthreads);
    nalloc_rb = (size_t*) malloc(sizeof(size_t)*nthreads);
    for (t = 0;t < nthreads;t++) {
        regionsbuffer[t] = (size_t*) malloc(sizeof(size_t) * 32);
        nalloc_rb[t] = 16;
    }
    nregions[0] = 1;
    regionsbuffer[0][0] = 0;
    regionsbuffer[0][1] = numel;
#pragma omp parallel
    {
        int cthread;
        int tthread;
        int stop;
        int hasnext;
        size_t begin;
        size_t end;
        size_t pivot;
        size_t cr;
        size_t startat;
        size_t ntodelegate;
        stop = ((*p_activethreads) == 0);
        cr = 0;
        while (!stop) {
            while (nregions[cthread] == 0 && !stop) {
                stop = ((*p_activethreads) == 0);
            }
            if (!stop) {
                begin = regionsbuffer[cthread][cr * 2];
                end = regionsbuffer[cthread][cr * 2 + 1];
                pivot = _vec_quicksortbyi64_makepivot(target, elt_size, values, begin, end);
                if (pivot < end) {
                    /* replace the current region with the region that's on the left of the pivot */
                    /* and add a new region with everything that's on the right of it *//*
                    if (nalloc_rb[cthread] == nregions[cthread]) {
                        nalloc_rb[cthread] *= 2;
                        regionsbuffer[cthread] = (size_t*)realloc(regionsbuffer[cthread], (sizeof(size_t)*nalloc_rb[cthread]*2));
                    }
                    regionsbuffer[cthread][nregions[cthread] * 2 + 1] = regionsbuffer[cthread][cr * 2 + 1];
                    regionsbuffer[cthread][cr * 2 + 1] = pivot;
                    regionsbuffer[cthread][nregions[cthread] * 2] = pivot;
#pragma omp critical
                    {
                        if ((*p_activethreads) < nthreads) {
                            /* delegate half of the remaining work to another thread *//*
                            *p_activethreads += 1;
                            for (tthread = 0;tthread < nthreads;tthread++) {
                                if (nregions[tthread] == 0)break;
                            }
                            startat = (nregions[cthread] + 1 - cr) / 2;
                            ntodelegate = nregions[cthread] - startat;
                            while (nalloc_rb[tthread] < ntodelegate) {
                                nalloc_rb[tthread] *= 2;
                                regionsbuffer[tthread] = (size_t*)realloc(regionsbuffer[tthread], (sizeof(size_t)*nalloc_rb[tthread] * 2));
                            }
                            memcpy(regionsbuffer[tthread], regionsbuffer[cthread] + startat, sizeof(size_t)*ntodelegate);
                            nregions[cthread] -= ntodelegate;
                            (*p_activethreads)++;
                            /* the following will trigger the target thread to stop waiting, so everything needs to be properly initialized at this point *//*
                            nregions[tthread] = ntodelegate;
                        }
                    }
                }
                else {
                    cr++;
                    if (cr >= nregions[cthread]) {
                        /* we have completed all the work on this segment *//*
                        cr = 0;
                        nregions[cthread] = 0;
#pragma omp atomic
                        (*p_activethreads)--;
                    }
                }
            }
        }
        free(regionsbuffer);
    }
    
}
/**/


uint64_t* vec_toranks(double* vec, size_t len) {
    /* this may not be the fastet way to get ranks - but it's easily implmented */
    uint64_t* order;
    uint64_t* result;
    double* tmpvec;
    uint64_t i;
    order = malloc(sizeof(uint64_t)*len);
    result = malloc(sizeof(uint64_t)*len);
    tmpvec = vec_cpy(vec, len);
    for (i = 0;i < (uint64_t)len;i++) {
        order[i] = i;
    }
    vec_sort_by((void*)order, sizeof(uint64_t), vec, len);
    for (i = 0;i < (uint64_t)len;i++) {
        result[order[i]] = i;
    }
    free(tmpvec);
    free(order);
    return result;
}
double vec_variance(double* X, size_t len) {
    double* tmpvec;
    double meanvalue;
    double squaresum;
    double result;
    tmpvec = vec_cpy(X, len);
    meanvalue = vec_avg(tmpvec, len);
    vec_dot(tmpvec, tmpvec, len);
    squaresum = vec_sum(tmpvec, len);
    result = squaresum / ((double)len) - meanvalue*meanvalue;
    free(tmpvec);
    return result;
}
double vec_variance_weighted(double* X, double* W, size_t len) {
    double* tmpvec;
    double meanvalue;
    double squaresum;
    double result;
    tmpvec = vec_cpy(X, len);
    meanvalue = vec_avg(tmpvec, len);
    vec_dot(tmpvec, tmpvec, len);
    vec_dot(tmpvec, W, len);
    squaresum = vec_sum(tmpvec, len);
    result = squaresum / ((double)len) - meanvalue*meanvalue;
    free(tmpvec);
    return result;
}
double vec_covariance(double* X, double* Y, size_t len) {
    double* tmpvec;
    double meanvalueX,meanvalueY;
    double squaresum;
    double result;
    meanvalueX = vec_avg(X, len);
    meanvalueY = vec_avg(Y, len);
    tmpvec = vec_cpy(X, len);
    vec_dot(tmpvec, Y, len);
    squaresum = vec_sum(tmpvec, len);
    result = squaresum / ((double)len) - meanvalueX*meanvalueY;
    free(tmpvec);
    return result;
}
double vec_pearsoncorr(double* X, double* Y, size_t len) {
    double* tmpvec;
    double meanvalueX, meanvalueY;
    double smvX, smvY; /* squared mean value of [] */
    double denomX, denomY;
    double squaresum;
    double numerator;
    meanvalueX = vec_avg(X, len);
    meanvalueY = vec_avg(Y, len);
    smvX = meanvalueX*meanvalueX;
    smvY = meanvalueY*meanvalueY;

    tmpvec = vec_cpy(X, len);
    vec_dot(tmpvec, Y, len);
    squaresum = vec_sum(tmpvec, len);
    numerator = squaresum/((double)len) - meanvalueX*meanvalueY;
   
    memcpy(tmpvec, X, sizeof(double)*len);
    vec_dot(tmpvec, tmpvec, len);
    denomX = vec_sum(tmpvec, len)/((double)len) - smvX;
    memcpy(tmpvec, Y, sizeof(double)*len);
    vec_dot(tmpvec, tmpvec, len);
    denomY = vec_sum(tmpvec, len)/((double)len) - smvY;

    free(tmpvec);
    return numerator / sqrt(denomX*denomY);
}
double vec_spearmancorr(double* X, double* Y, size_t len) {
    uint64_t* ranksX;
    uint64_t* ranksY;
    double* dranksX;
    double* dranksY;
    double result;

    ranksX = vec_toranks(X, len);
    ranksY = vec_toranks(Y, len);
    dranksX = veci64_to_vec(ranksX, len);
    free(ranksX);
    dranksY = veci64_to_vec(ranksY, len);
    free(ranksY);

    result = vec_pearsoncorr(dranksX, dranksY, len);

    free(dranksX);
    free(dranksY);

    return result;
}
double vec_manhattandist(double* A, double* B, size_t len) {
    double result;
    double delta;
    size_t i;
    result = 0.0;
    for (i = 0;i < len;i++) {
        delta = A[i] - B[i];
        if (delta > 0) result += delta;
        else result -= delta;
    }
    return result;
}
double vec_euclidiandist(double* A, double* B, size_t len) {
    double result;
    double delta;
    size_t i;
    result = 0.0;
    for (i = 0;i < len;i++) {
        delta = A[i] - B[i];
        result += delta * delta;
    }
    return sqrt(result);
}
