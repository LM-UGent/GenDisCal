#include <math.h>
#include <string.h>
#include "datatable.h"
#include "datamap.h"

struct datatable_t {
    double** values;
    size_t w;
    size_t h;
    DM64_t* rowids;
    DM64_t* colids;
    char** rownames;
    char** colnames;
    double fill;
};

datatable_t* datatable_alloc(size_t ncol, size_t nrow, double autofill) {
    datatable_t* result;
    size_t i,j;
    result = (datatable_t*) calloc(1, sizeof(datatable_t));
    result->w = ncol;
    result->h = nrow;
    result->values = (double**)calloc(nrow, sizeof(double*));
    for (i = 0;i < nrow;i++) {
        result->values[i] = (double*)malloc(ncol*sizeof(double));
        for (j = 0;j < ncol;j++) {
            result->values[i][j] = autofill;
        }
    }
    result->rowids = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
    result->colids = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
    result->rownames = (char**)calloc(nrow, sizeof(char*));
    result->colnames = (char**)calloc(ncol, sizeof(char*));
    result->fill = autofill;
    return result;
}
void datatable_free(datatable_t* target) {
    size_t i;
    for (i = 0;i < target->h;i++) {
        free(target->values[i]);
        target->values[i] = NULL;
        if (target->rownames[i])free(target->rownames[i]);
        target->rownames[i] = NULL;
    }
    free(target->values);
    free(target->rownames);
    for (i = 0;i < target->w;i++) {
        if (target->colnames[i])free(target->colnames[i]);
        target->colnames[i] = NULL;
    }
    free_DM64(target->colids);
    free_DM64(target->rowids);
    target->colids = NULL;
    target->rowids = NULL;
}

size_t datatable_ncols(datatable_t* target) {
    return target->w;
}
size_t datatable_nrows(datatable_t* target){
    return target->h;
}

void _datatable_addcolsupto(datatable_t* target, size_t col) {
    size_t i;
    size_t j;
    if (col >= target->w) {
        for (i = 0;i < target->h;i++) {
            target->values[i] = (double*)realloc(target->values[i], sizeof(double)*(col + 1));
            for (j = target->w; j <= col;j++) {
                target->values[i][j] = target->fill;
            }
        }
        target->colnames = realloc(target->colnames, sizeof(char*)*(col + 1));
        for (j = target->w; j <= col;j++) {
            target->colnames[j] = NULL;
        }
        target->w = col + 1;
    }
}
void _datatable_addrowsupto(datatable_t* target, size_t row) {
    size_t i;
    size_t j;
    if(row >= target->h) {
        target->values = (double**)realloc(target->values, sizeof(double*)*(row + 1));
        target->rownames = (char**)realloc(target->rownames, sizeof(char*)*(row + 1));
        for (i = target->h;i <= row;i++) {
            target->values[i] = (double*)malloc(sizeof(double)*target->w);
            target->rownames[i] = NULL;
            for (j = 0; j < target->w;j++) {
                target->values[i][j] = target->fill;
            }
        }
        target->h = row + 1;
    }
}

errcode_t datatable_set(datatable_t* target, size_t col, size_t row, double value) {
    _datatable_addcolsupto(target, col);
    _datatable_addrowsupto(target, row);
    target->values[row][col] = value;
    return 0;
}
double datatable_get(datatable_t* target, size_t col, size_t row, int* nullflag){
    if (col >= target->w || row >= target->h) {
        *nullflag = 1;
        return target->fill;
    }
    return target->values[row][col];
}

void datatable_setcolname(datatable_t* target, size_t col, const char* name) {
    size_t namelen;
    _datatable_addcolsupto(target, col);
    if (target->colnames[col]) free(target->colnames[col]);
    namelen = strlen(name) + 1;
    target->colnames[col] = malloc(namelen);
    memcpy(target->colnames[col], name, namelen);
    DM64_append(target->colids, name, (int)namelen - 1, (int64_t)col);
}
size_t datatable_getcolid(datatable_t* target, const char* name){
    int nullflag;
    int64_t result;
    result = DM64_get(target->colids, name, (int) strlen(name), &nullflag);
    if (nullflag)
        result = target->w;
    return result;
}
char* datatable_getcolname(datatable_t* target, size_t col) {
    if (col >= target->w)return NULL;
    return target->colnames[col];
}

void datatable_setrowname(datatable_t* target, size_t row, const char* name) {
    size_t namelen;
    _datatable_addrowsupto(target, row);
    if (target->rownames[row]) free(target->rownames[row]);
    namelen = strlen(name) + 1;
    target->rownames[row] = malloc(namelen);
    memcpy(target->rownames[row], name, namelen);
    DM64_append(target->rowids, name, (int)namelen - 1, (int64_t)row);
}
ssize_t datatable_getrowid(datatable_t* target, const char* name) {
    int nullflag;
    int64_t result;
    result = DM64_get(target->rowids, name, (int) strlen(name), &nullflag);
    if (nullflag)
        result = target->h;
    return result;
}
char* datatable_getrowname(datatable_t* target, size_t row) {
    if (row >= target->h)return NULL;
    return target->rownames[row];
}

void datatable_write(datatable_t* dt, PF_t* target, const char* sep, uint32_t flags) {
    size_t i, j;
    if (flags & DATATABLE_WRITECOLNAMES) {
        if (flags & DATATABLE_WRITEROWNAMES) {
            PFprintf(target, "%s", sep);
        }
        if (dt->w > 0) {
            for (j = 0;j < dt->w - 1;j++) {
                if(dt->colnames[j])
                    PFprintf(target, "%s", dt->colnames[j]);
                PFprintf(target, "%s", sep);
            }
            if (dt->colnames[j])
                PFprintf(target, "%s", dt->colnames[j]);
        }
        PFprintf(target, "\n");
    }
    for (i = 0;i < dt->h;i++) {
        if (flags & DATATABLE_WRITEROWNAMES) {
            if (dt->rownames[i])
                PFprintf(target, "%s", dt->rownames[i]);
            PFprintf(target, "%s", sep);
        }
        for (j = 0;j < dt->w - 1;j++) {
            for (j = 0;j < dt->w - 1;j++) {
                PFprintf(target, "%f%s", dt->values[i][j], sep);
            }
            if (dt->colnames[j])
                PFprintf(target, "%f", dt->values[i][j]);
        }
        PFprintf(target, "\n");
    }
}