#ifndef DEF_HASH_TABLE
#define DEF_HASH_TABLE

#include "osportstd.h"
#include "hash.h"

typedef struct hash_table_32_t ht32_t;
typedef struct hash_table_64_t ht64_t;

ht32_t* ht32_alloc();
ht32_t* ht32_alloc_size(size_t numel);
void ht32_free(ht32_t* target);
int32_t ht32_get(ht32_t* target, char* key, size_t keylen, int* nullflag);
void ht32_set(ht32_t* target, char* key, size_t keylen, int32_t value, int* nullflag);
void ht32_inc(ht32_t* target, char* key, size_t keylen, int32_t by, int* nullflag);
size_t ht32_astables(ht32_t * target, char*** p_listofnames, size_t** p_listoflens, int32_t** p_values);


ht64_t* ht64_alloc();
ht64_t* ht64_alloc_size(size_t numel);
void ht64_free(ht64_t* target);
int64_t ht64_get(ht64_t* target, char* key, size_t keylen, int* nullflag);
void ht64_set(ht64_t* target, char* key, size_t keylen, int64_t value, int* nullflag);
void ht64_inc(ht64_t* target, char* key, size_t keylen, int64_t by, int* nullflag);
size_t ht64_astables(ht64_t * target, char*** p_listofnames, size_t** p_listoflens, int64_t** p_values);
void ht64_save(ht64_t* table, PF_t* file);
ht64_t* ht64_load(PF_t* file);
#endif
