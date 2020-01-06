#ifndef DEF_DATATABLE
#define DEF_DATATABLE

#include "osportstd.h"

typedef struct datatable_t datatable_t;

typedef struct dendrogram_t dendrogram_t;

datatable_t* datatable_alloc(size_t ncol, size_t nrow, double autofill);
void datatable_free(datatable_t* target);

size_t datatable_ncols(datatable_t* target);
size_t datatable_nrows(datatable_t* target);

errcode_t datatable_set(datatable_t* target, size_t col, size_t row, double value);
double datatable_get(datatable_t* target, size_t col, size_t row, int* nullflag);

void datatable_addtables(datatable_t* target, datatable_t** toadd, size_t numtables);

void datatable_setcolname(datatable_t* target, size_t col, const char* name);
size_t datatable_getcolid(datatable_t* target, const char* name);
char* datatable_getcolname(datatable_t* target, size_t col);

void datatable_setrowname(datatable_t* target, size_t col, const char* name);
ssize_t datatable_getrowid(datatable_t* target, const char* name);
char* datatable_getrowname(datatable_t* target, size_t col);

/*
dmrclust stands for distance matrix rough clustering, and will only work on distance matrices.
Also, it is a rough sorting algorithm with complexity k*n*log(n) which focuses on clustering
all closely related vectors together, but not ensureing that proper ordering is observed
within clusters. Another algorithm should be used on the clusters to achive this.
*/
size_t datatable_dmrclust(datatable_t* target, double defaultcut);
size_t datatable_dmwardclust(datatable_t* target, double defaultcut, PF_t* output_tree);
size_t datatable_dmsingleclust(datatable_t* target, PF_t* output_tree);
size_t datatable_dmsortclust(datatable_t* target, double defaultcut);
#define DATATABLE_WRITECOLNAMES     1
#define DATATABLE_WRITEROWNAMES     2
void datatable_write(datatable_t* dt, PF_t* target, const char* sep, uint32_t flags);


#endif

