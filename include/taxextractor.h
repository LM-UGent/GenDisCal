#ifndef DEF_TAXEXTRACTOR
#define DEF_TAXEXTRACTOR

#include <stdlib.h>
#include <stdint.h>

typedef struct taxextractor_t taxextractor_t;

taxextractor_t* taxextractor_alloc();
void taxextractor_free(taxextractor_t* target);
char** _taxextractor_PCOF_S_to_str_array(const char* string, size_t* count);
void taxextractor_add_PCOF_S(taxextractor_t* extractor, const char* string);
int taxextractor_translate_PCOF_S(taxextractor_t* extractor, const char* string, int64_t* idarray, size_t idarraysize);
int taxextractor_add(taxextractor_t* extractor, const char* method, const char* string);
int taxextractor_translate(taxextractor_t* extractor, const char* method, const char* string, int64_t* idarray, size_t idarraysize);
void taxextractor_sort(taxextractor_t* extractor);
#endif

