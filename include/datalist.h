#ifndef DEF_DATALIST
#define DEF_DATALIST

#include "stdint.h"

int64_t ___added_values___(int64_t delta);

typedef struct doublelink64_t {
    int64_t value;
    struct doublelink64_t* next;
    struct doublelink64_t* prev;
} dlink64_t;

dlink64_t* dlink_alloc(int64_t value);
dlink64_t* dlink_insert_after(dlink64_t* at, int64_t value);
dlink64_t* dlink_insert_before(dlink64_t* at, int64_t value);
dlink64_t* dlink_append(dlink64_t* list, int64_t value);
dlink64_t* dlink_remove(dlink64_t* target);
void dlink_freeall(dlink64_t* target);


typedef struct heaplink64_t {
    int64_t value;
    struct heaplink64_t* parent;
    struct heaplink64_t* lchild;
    struct heaplink64_t* rchild;
} hlink64_t;

hlink64_t* hlink_insert(hlink64_t* root, int64_t value);
hlink64_t* hlink_remove(hlink64_t* target);
void hlink_freeall(hlink64_t* target);


#endif
