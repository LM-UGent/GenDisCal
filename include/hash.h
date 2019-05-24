#ifndef DEF_HASH
#define DEF_HASH

#include <stdint.h>
#include <stddef.h>

typedef struct hash_descriptor_t hashdesc_t;

hashdesc_t* hashdesc_alloc();
void hashdesc_free(hashdesc_t* target);
int hashdesc_init_basic(hashdesc_t* target, size_t maxelt);
int hashdesc_init_fingerprint64(hashdesc_t* target);
size_t hash_index(void* data, size_t datasize, hashdesc_t* desc);



#endif
