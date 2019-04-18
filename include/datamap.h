#ifndef DEF_DATAMAP
#define DEF_DATAMAP

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct _datamap32 DM32_t;
typedef struct _datamap64 DM64_t;

DM32_t* new_DM32(int algorithm, int32_t default_value);
DM64_t* new_DM64(int algorithm, int64_t default_value);

void free_DM32(DM32_t* target);
void free_DM64(DM64_t* target);

int DM32_write(DM32_t* target, FILE* f);
int DM64_write(DM64_t* target, FILE* f);

int DM32_read(DM32_t* target, FILE* f);
int DM64_read(DM64_t* target, FILE* f);

#define DM_ALGORITHM_BASICSORTEDLIST   0   /* a basic and inefficient algorithm, but simple to implement */
/*
appending is preferable when adding many elements to the datamap as it does not change the internal structure
the DMX_append functions do just that, but this leaves the datamap unusable. The DMX_sort methods can put the
datamap back in the proper order, and are called automatically at the next get operation.
*/
void DM32_append(DM32_t* datamap, const char* key, int keylen, int32_t value);
void DM64_append(DM64_t* datamap, const char* key, int keylen, int64_t value);

void DM32_sort(DM32_t* datamap);
void DM64_sort(DM64_t* datamap);

void DM32_assign(DM32_t* datamap, const char* key, int keylen, int32_t value);
void DM64_assign(DM64_t* datamap, const char* key, int keylen, int64_t value);

int32_t DM32_get(DM32_t* datamap, const char* key, int keylen, int* NULLflag);
int64_t DM64_get(DM64_t* datamap, const char* key, int keylen, int* NULLflag);

void DM32_printf(DM32_t* datamap);
void DM64_printf(DM64_t* datamap);

size_t DM32_numel(DM32_t* datamap);
size_t DM64_numel(DM64_t* datamap);

void DM32_ordertovalues(DM32_t* datamap);
void DM64_ordertovalues(DM64_t* datamap);

#endif // DEF_DATAMAP
