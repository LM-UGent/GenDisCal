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

#include "hash.h"
#include <stdlib.h>


/* TOTO: when a new hash is implemented, add it here*/
typedef enum hash_type { BASIC_HASH = 0 } hash_type;

struct hash_descriptor_t {
    hash_type id;
    size_t maxavlue;
    uint64_t v1;
    uint64_t v2;
    uint64_t v3;
};

hashdesc_t* hashdesc_alloc() {
    return (hashdesc_t*)calloc(1, sizeof(hashdesc_t));
}
void hashdesc_free(hashdesc_t* target) {
    if(target) free(target);
}
int hashdesc_init_basic(hashdesc_t* target, size_t maxelt) {
    target->id = BASIC_HASH;
    target->maxavlue = maxelt;
    if (maxelt < 0x0fffffff) {
        /* should be sufficient for most use cases */
        target->v1 = 48271;
        target->v2 = 0x7fffffff;
        target->v3 = 0;
    }
    else if(maxelt < 0x1000000000000){
        /* Extend the limit to 2^48 (~281T elements), which should be sufficient in all reasonable cases */
        target->v1 = 25214903917;
        target->v2 = 0x1000000000000;
        target->v3 = 11;
    }
    else {
        /* The descriptor could not be initialized due to unreasonable parameters */
        return 1;
    }
    return 0;
}
int hashdesc_init_fingerprint64(hashdesc_t* target) {
    return hashdesc_init_basic(target, 0xFFFFFFFFFFFF);
}
size_t _hash_index_basic(void* data, size_t datasize, hashdesc_t* desc) {
    /*
    The "basic" hash works essentially in the same way as the minstd random
    number generator. This should be sufficient for most applications.
    */
    size_t i;
    uint64_t h;
    uint64_t factor;
    char* x;
    h = datasize*desc->v3;
    x = (char*)data;
    factor = desc->v1;
    for (i = 0; i < datasize; i++)
        h = ((h*factor) + x[i] + desc->v3) % desc->v2;
    h = ((h*factor) + 1) % desc->v2;
    return h%desc->maxavlue;
}

size_t hash_index(void* data, size_t datasize, hashdesc_t* desc) {
    if (!data)return -1;
    if (desc->id == BASIC_HASH) return _hash_index_basic(data, datasize, desc);
    return -1;
}

void hashdesc_save(hashdesc_t* desc, PF_t* file) {
    PFputint32(file, (int32_t)(desc->id));
    PFputint64(file, (int64_t)(desc->maxavlue));
    PFputint64(file, desc->v1);
    PFputint64(file, desc->v2);
    PFputint64(file, desc->v3);
}
hashdesc_t* hashdesc_load(PF_t* file) {
    hashdesc_t* result;
    result = hashdesc_alloc();
    result->id = (hash_type)PFgetint32(file);
    result->maxavlue = (size_t)PFgetint64(file);
    result->v1 = PFgetint64(file);
    result->v2 = PFgetint64(file);
    result->v3 = PFgetint64(file);
    return result;
}
