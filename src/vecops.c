#include "vecops.h"

double* vec_cpy(double* source, size_t len) {
    double* result;
    size_t i = 0;
    result = malloc(sizeof(double)*len);
    for (i = 0;i < len;i++) result[i] = source[i];
    return result;
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
    //not implemented
    return 0.0;
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

#define arrswapids(arr,index1,index2,tmpvar)     tmpvar = arr[index1]; arr[index1] = arr[index2]; arr[index2] = tmpvar
void _vec_quicksort(double* target, size_t len) {
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
void _vec_heapify(double* target, size_t len, size_t node) {
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
void _vec_sort__wdepth(double* target, size_t len, size_t depth, size_t maxdepth) {
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
            /* doing this prevents an infinite loop */
        }
        _vec_quicksort(target, pivot);
        _vec_quicksort(target, len - pivot);
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
