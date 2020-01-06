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

#include "hash_table.h"
#include <stdlib.h>
#include <string.h>

typedef char* keyptr_t;

static inline size_t _htkeydata_length(keyptr_t key) {
    return (size_t)((size_t*)(key+1))[0];
}
static inline int _htkeydata_type(keyptr_t key) {
    return (int)key[0];
}
static inline char* _htkeydata_key(keyptr_t key) {
    return key + sizeof(size_t)+1;
}
static inline keyptr_t _htkeydata_nextkey(keyptr_t key, size_t extrasize) {
    return ((keyptr_t*)(key+1+sizeof(size_t)+_htkeydata_length(key)+extrasize))[0];
}
static inline void _htkeydata_setnext(keyptr_t prevkey, keyptr_t nextkey, size_t extrasize) {
    prevkey[0] = (char)2;
    nextkey[0] = (char)1;
    ((keyptr_t*)(prevkey + 1 + sizeof(size_t) + _htkeydata_length(prevkey) + extrasize))[0] = nextkey;
}
static inline int32_t _htkeydata_value32(keyptr_t key) {
    return ((int32_t*)(key + _htkeydata_length(key) + 1 + sizeof(size_t)))[0];
}
static inline int64_t _htkeydata_value64(keyptr_t key) {
    return ((int64_t*)(key + _htkeydata_length(key) + 1 + sizeof(size_t)))[0];
}
static inline void* _htkeydata_valuep(keyptr_t key) {
    return (void*)(key + _htkeydata_length(key) + 1 + sizeof(size_t));
}
static inline void _htkeydata_setvalue32(keyptr_t key, int32_t value) {
    int32_t* valuep;
    size_t oldkeylen;
    oldkeylen = _htkeydata_length(key);
    valuep = (int32_t*)(key + 1 + oldkeylen + sizeof(size_t));
    *valuep = value;
}
static inline void _htkeydata_setvalue64(keyptr_t key, int64_t value) {
    int64_t* valuep;
    size_t oldkeylen;
    oldkeylen = _htkeydata_length(key);
    valuep = (int64_t*)(key + 1 + oldkeylen + sizeof(size_t));
    *valuep = value;
}
static inline keyptr_t _htkeydata_create(char* keyvalue, size_t len, size_t extrasize) {
    char* result;
    size_t* lenp;
    char* keyp;
    result = malloc(len + sizeof(size_t) + 1 + extrasize);
    result[0] = 1;
    lenp = (size_t*)(result + 1);
    *lenp = len;
    keyp = result + 1 + sizeof(size_t);
    memcpy(keyp, keyvalue, len);
    return result;
}
static keyptr_t _htkeydata_create32(char* keyvalue, size_t len, int32_t value) {
    char* result;
    int32_t* pvalue;
    result = _htkeydata_create(keyvalue, len, sizeof(int32_t));
    pvalue = (int32_t*)(result + 1 + sizeof(size_t) + len);
    *pvalue = value;
    return result;
}
static keyptr_t _htkeydata_create64(char* keyvalue, size_t len, int64_t value) {
    char* result;
    int64_t* pvalue;
    result = _htkeydata_create(keyvalue, len, sizeof(int64_t));
    pvalue = (int64_t*)(result + 1 + sizeof(size_t) + len);
    *pvalue = value;
    return result;
}
static keyptr_t _htkeydata_getkeyptr(keyptr_t* keyaddr, char* keyvalue, size_t len, size_t* depth, size_t extrasize, int* nullflag) {
    int keytype;
    size_t oldkeylen;
    char* oldkeyname;
    char* result;
    size_t oldend;
    size_t tmpdepth;
    keyptr_t* newkey;
    keyptr_t key;
    keyptr_t* parentkey;
    *depth = 0;
    parentkey = keyaddr;
    if (!keyaddr) return NULL;
    key = *keyaddr;
    if (!key)return NULL;

    keytype = _htkeydata_type(key);
    oldkeylen = _htkeydata_length(key);
    oldkeyname = _htkeydata_key(key);
    result = NULL;
    tmpdepth = 0;
    *nullflag = 0;
    while (keytype == 2 && result==NULL) {
        if (oldkeylen == len && memcmp(oldkeyname, keyvalue, len)==0)
            result = key;
        else {
            /* get a pointer to the adress value stored at the end of the previous key */
            parentkey = (keyptr_t*)(key+extrasize+1+oldkeylen+sizeof(size_t));
            key = _htkeydata_nextkey(key, extrasize);
            keytype = _htkeydata_type(key);
            oldkeylen = _htkeydata_length(key);
            oldkeyname = _htkeydata_key(key);
        }
        tmpdepth++;
    }
    if (oldkeylen == len && memcmp(oldkeyname, keyvalue, len) == 0)
        result = key;
    else if (keytype == 1 && result==NULL) {
        *nullflag = 1;
        oldend = 1 + oldkeylen + sizeof(size_t) + extrasize;
        tmpdepth++;
        key = realloc(key, oldend + sizeof(char*));
        *parentkey = key;
        key[0] = 2;
        newkey = (keyptr_t*)(key + oldend);
        *newkey = _htkeydata_create(keyvalue, len, extrasize);
        result = *newkey;
    }
    *depth = tmpdepth;
    return result;
}
static keyptr_t _htkeydata_set32(keyptr_t key, char* keyvalue, size_t len, size_t* depth, int32_t value, int* nullflag) {
    char* target;
    target = _htkeydata_getkeyptr(&key, keyvalue, len, depth, sizeof(int32_t), nullflag);
    _htkeydata_setvalue32(target, value);
    return key;
}
static keyptr_t _htkeydata_set64(keyptr_t key, char* keyvalue, size_t len, size_t* depth, int64_t value, int* nullflag) {
    char* target;
    target = _htkeydata_getkeyptr(&key, keyvalue, len, depth, sizeof(int64_t), nullflag);
    _htkeydata_setvalue64(target, value);
    return key;
}
static keyptr_t _htkeydata_rm(keyptr_t key, char* keyvalue, size_t len, size_t extrasize) {
    int keytype;
    size_t oldkeylen;
    char* oldkeyname;
    char* result;
    char* targetkey;
    char* parentkey;
    char** nextkeyp;
    char* childkey;
    if (!key)return NULL;

    keytype = _htkeydata_type(key);
    oldkeylen = _htkeydata_length(key);
    oldkeyname = _htkeydata_key(key);
    result = NULL;
    parentkey = NULL;
    childkey = NULL;
    targetkey = key;
    /* find a matching key */
    while (keytype == 2 && result == NULL) {
        if (oldkeylen == len && memcmp(oldkeyname, keyvalue, len)==0)
            result = targetkey;
        else {
            parentkey = targetkey;
            targetkey = _htkeydata_nextkey(targetkey, extrasize);
            keytype = _htkeydata_type(targetkey);
            oldkeylen = _htkeydata_length(targetkey);
            oldkeyname = _htkeydata_key(targetkey);
        }
    }
    if (oldkeylen == len && memcmp(oldkeyname, keyvalue, len) == 0)
        result = targetkey;
    if (result){
        /* if a matching key was found, remove it*/
        if (keytype == 2) {
            childkey = _htkeydata_nextkey(result, extrasize);
        }
        if (!parentkey) {
            parentkey = childkey;
            key = childkey;
        }
        else {
            /* make the parentkey point to the child key*/
            oldkeylen = _htkeydata_length(parentkey);
            nextkeyp = (char**)(parentkey + 1 + sizeof(size_t) + extrasize + oldkeylen);
            *nextkeyp = childkey;
            if (!childkey) {
                /* if parentkey is the now the last node, change its type to reflect this*/
                parentkey[0] = 1;
            }
        }
    }
    return key;
}
static keyptr_t _htkeydata_rm32(keyptr_t key, char* keyvalue, size_t len) {
    return _htkeydata_rm(key, keyvalue, len, sizeof(int32_t));
}
static keyptr_t _htkeydata_rm64(keyptr_t key, char* keyvalue, size_t len) {
    return _htkeydata_rm(key, keyvalue, len, sizeof(int64_t));
}
static void _htkeydata_free(keyptr_t key, size_t extrasize) {
    int keytype;
    size_t oldkeylen;
    char* oldkeyname;
    char* nextkey;
    if (!key)return;

    keytype = _htkeydata_type(key);
    oldkeylen = _htkeydata_length(key);
    oldkeyname = _htkeydata_key(key);
    /* find a matching key */
    while (keytype == 2) {
        nextkey = _htkeydata_nextkey(key, extrasize);
        free(key);
        key = nextkey;
        keytype = _htkeydata_type(key);
        oldkeylen = _htkeydata_length(key);
        oldkeyname = _htkeydata_key(key);
    }
    free(key);
}
static void _htkeydata_free32(keyptr_t key) {
    _htkeydata_free(key, sizeof(int32_t));
}
static void _htkeydata_free64(keyptr_t key) {
    _htkeydata_free(key, sizeof(int64_t));
}

struct hash_table_32_t {
    hashdesc_t* hashdesc;
    size_t nalloc;
    keyptr_t* keys;
    size_t nset;
};

ht32_t* ht32_alloc() {
    ht32_t* result;
    result = malloc(sizeof(ht32_t));

    result->nalloc = 16;
    result->keys = (keyptr_t*)calloc(16, sizeof(keyptr_t));
    result->hashdesc = hashdesc_alloc();
    result->nset = 0;
    hashdesc_init_basic(result->hashdesc, 16);
    return result;
}
ht32_t* ht32_alloc_size(size_t numel)
{
    ht32_t* result;
    result = malloc(sizeof(ht32_t));

    result->nalloc = numel;
    result->keys = (keyptr_t*)calloc(numel, sizeof(keyptr_t));
    result->hashdesc = hashdesc_alloc();
    result->nset = 0;
    hashdesc_init_basic(result->hashdesc, numel);
    return result;
}
void ht32_free(ht32_t* target) {
    size_t i;
    if (target->nalloc) {
        for (i = 0;i < target->nalloc;i++) {
            if (target->keys[i])
                _htkeydata_free32(target->keys[i]);
        }
        free(target->keys);
        target->nalloc = 0;
    }
    hashdesc_free(target->hashdesc);
    free(target);
}
int32_t ht32_get(ht32_t* target, char* key, size_t keylen, int* nullflag) {
    size_t hash_val;
    keyptr_t keyptr;
    int recorded_key_type;
    size_t recorded_key_length;
    int32_t result;
    int found;
    *nullflag = 1;
    if (!target)return 0;
    result = 0;
    hash_val = hash_index(key, keylen, target->hashdesc);
    if (target->keys[hash_val]) {
        keyptr = target->keys[hash_val];
        recorded_key_type = _htkeydata_type(keyptr);
        recorded_key_length = _htkeydata_length(keyptr);
        if (recorded_key_length == keylen && memcmp(key, _htkeydata_key(keyptr), keylen) == 0)
            found = 1;
        else
            found = 0;
        while (recorded_key_type == 2 && !found) {
            keyptr = _htkeydata_nextkey(keyptr, sizeof(int32_t));
            recorded_key_type = _htkeydata_type(keyptr);
            recorded_key_length = _htkeydata_length(keyptr);
            if (recorded_key_length == keylen && memcmp(key, _htkeydata_key(keyptr), keylen) == 0)
                found = 1;
        }
        if (found) {
            *nullflag = 0;
            result = _htkeydata_value32(keyptr);
        }
    }
    return result;
}
static void _ht32_rehash(ht32_t* target, size_t newnumel) {
    hashdesc_t* newhdesc;
    size_t i;
    keyptr_t* newkeys;
    keyptr_t tmpkey;
    keyptr_t tmpkey2;
    size_t keylen;
    int32_t keyvalue;
    size_t newindex;
    size_t depth;
    char* keyname;
    int nullflag;
    newhdesc = target->hashdesc;
    hashdesc_init_basic(newhdesc, newnumel);
    newkeys = (keyptr_t*)calloc(newnumel, sizeof(keyptr_t));
    for (i = 0;i < target->nalloc;i++) {
        if (target->keys[i]) {
            tmpkey = target->keys[i];
            while (_htkeydata_type(tmpkey) == 2)
            {
                keylen = _htkeydata_length(tmpkey);
                keyname = _htkeydata_key(tmpkey);
                keyvalue = _htkeydata_value32(tmpkey);
                newindex = hash_index((void*)keyname, keylen, newhdesc);
                if (newkeys[newindex] == NULL) {
                    newkeys[newindex] = tmpkey;
                    ((char*)tmpkey)[0] = 1; /* change the type to be an end node*/
                    tmpkey = _htkeydata_nextkey(tmpkey, sizeof(int32_t));
                }
                else {
                    newkeys[newindex] = _htkeydata_set32(newkeys[newindex], keyname, keylen, &depth, keyvalue, &nullflag);
                    tmpkey2 = _htkeydata_nextkey(tmpkey, sizeof(int32_t));
                    _htkeydata_rm(tmpkey, keyname, keylen, sizeof(int32_t));
                    tmpkey = tmpkey2;
                }
            }
            keylen = _htkeydata_length(tmpkey);
            keyname = _htkeydata_key(tmpkey);
            keyvalue = _htkeydata_value32(tmpkey);
            newindex = hash_index((void*)keyname, keylen, newhdesc);
            if (newkeys[newindex] == NULL) {
                newkeys[newindex] = tmpkey;
                ((char*)tmpkey)[0] = 1; /* change the type to be an end node*/
            }
            else {
                newkeys[newindex] = _htkeydata_set32(newkeys[newindex], keyname, keylen, &depth, keyvalue, &nullflag);
                tmpkey2 = _htkeydata_nextkey(tmpkey, sizeof(int32_t));
                _htkeydata_rm(tmpkey, keyname, keylen, sizeof(int32_t));
                tmpkey = tmpkey2;
            }
        }
    }
    free(target->keys);
    target->keys = newkeys;
    target->nalloc = newnumel;
}
void ht32_set(ht32_t* target, char* key, size_t keylen, int32_t value, int *nullflag) {
    size_t hash_val;
    keyptr_t keyptr;
    size_t depth;
    *nullflag = 1;
    if (!target)return;
    hash_val = hash_index(key, keylen, target->hashdesc);
    if (target->keys[hash_val]) {
        keyptr = target->keys[hash_val];
        target->keys[hash_val] = _htkeydata_set32(keyptr, key, keylen, &depth, value, nullflag);
        if (*nullflag) {
            target->nset++;
            if (target->nset > target->nalloc)
                _ht32_rehash(target, target->nalloc * 2);
        }
    }
    else {
        target->keys[hash_val] = _htkeydata_create32(key, keylen, value);
        *nullflag = 1;
        target->nset++;
    }
}
void ht32_inc(ht32_t * target, char * key, size_t keylen, int32_t by, int * nullflag)
{
    size_t hash_val;
    keyptr_t keyptr;
    size_t depth;
    int32_t oldval;
    *nullflag = 1;
    if (!target)return;
    hash_val = hash_index(key, keylen, target->hashdesc);
    if (target->keys[hash_val]) {
        keyptr = target->keys[hash_val];
        oldval = _htkeydata_value32(keyptr);
        target->keys[hash_val] = _htkeydata_set32(keyptr, key, keylen, &depth, oldval + by, nullflag);
        if (*nullflag) {
            target->nset++;
            if (target->nset > target->nalloc)
                _ht32_rehash(target, target->nalloc * 2);
        }
    }
    else {
        target->keys[hash_val] = _htkeydata_create32(key, keylen, by);
        *nullflag = 1;
        target->nset++;
    }
}
size_t ht32_astables(ht32_t * target, char*** p_listofnames,size_t** p_listoflens, int32_t** p_values)
{
    char** result;
    size_t* resultlen;
    int32_t* resultvalues;
    size_t i,j;
    size_t len;
    keyptr_t curkey;

    if (p_values)
        resultvalues = (int32_t*)malloc(target->nset * sizeof(int32_t));
    else
        resultvalues = NULL;
    if (p_listofnames)
        result = (char**)malloc(target->nset * sizeof(char*));
    else
        result = NULL;
    if (p_listoflens)
        resultlen = (size_t*)malloc(target->nset * sizeof(size_t));
    else
        resultlen = NULL;
    j = 0;
    curkey = target->keys[0];
    for (i = 0;i < target->nset;i++) {
        if (target->keys[j] == NULL) {
            while (target->keys[j] == NULL)
                j++;
            curkey = target->keys[j];
        }
        len = _htkeydata_length(curkey);
        if (result) {
            result[i] = malloc((len + 1) * sizeof(char));
            memcpy(result[i], _htkeydata_key(curkey), len);
            result[i][len] = 0;
        }
        if (resultlen)
            resultlen[i] = len;
        if (resultvalues)
            resultvalues[i] = _htkeydata_value32(curkey);
        if (_htkeydata_type(curkey) == 2)
            curkey = _htkeydata_nextkey(curkey, sizeof(int32_t));
        else {
            j++;
            curkey = target->keys[j];
        }
    }
    *p_listofnames = result;
    *p_values = resultvalues;
    *p_listoflens = resultlen;
    return i;
}
void ht32_save(ht32_t* table, PF_t* file) {
    size_t i, j;
    size_t len;
    keyptr_t curkey;
    j = 0;
    curkey = table->keys[0];
    PFputint64(file, (uint64_t)(table->nalloc));
    PFputint64(file, (uint64_t)(table->nset));
    hashdesc_save(table->hashdesc, file);
    for (i = 0;i < table->nset;i++) {
        if (table->keys[j] == NULL) {
            while (table->keys[j] == NULL)
                j++;
            curkey = table->keys[j];
        }
        PFputint64(file, (int64_t)j);
        len = _htkeydata_length(curkey);
        PFputint64(file, (int64_t)len);
        PFwrite(_htkeydata_key(curkey), 1, len, file);
        PFputint32(file, _htkeydata_value32(curkey));
        if (_htkeydata_type(curkey) == 2)
            curkey = _htkeydata_nextkey(curkey, sizeof(int32_t));
        else {
            j++;
            curkey = table->keys[j];
        }
    }
}
ht32_t* ht32_load(PF_t* file) {
    ht32_t* result;
    size_t i, j;
    size_t len;
    char* key;
    keyptr_t keyptr;
    size_t depth;
    int nf;
    int32_t value;
    result = malloc(sizeof(ht32_t));
    result->nalloc = (size_t)PFgetint64(file);
    result->nset = (size_t)PFgetint64(file);
    result->hashdesc = hashdesc_load(file);
    result->keys = (keyptr_t*)calloc(result->nalloc, sizeof(keyptr_t));
    for (i = 0;i < result->nset;i++) {
        j = (size_t)PFgetint64(file);
        len = (size_t)PFgetint64(file);
        value = PFgetint32(file);
        key = malloc(len);
        PFread(key, 1, len, file);
        if (result->keys[j]) {
            keyptr = result->keys[j];
            result->keys[j] = _htkeydata_set32(keyptr, key, len, &depth, value, &nf);
        }
        else {
            result->keys[j] = _htkeydata_create32(key, len, value);
        }
    }
    return result;
}



struct hash_table_64_t {
    hashdesc_t* hashdesc;
    size_t nalloc;
    keyptr_t* keys;
    size_t nset;
};

ht64_t* ht64_alloc() {
    ht64_t* result;
    result = malloc(sizeof(ht64_t));

    result->nalloc = 16;
    result->keys = (keyptr_t*)calloc(16, sizeof(keyptr_t));
    result->hashdesc = hashdesc_alloc();
    result->nset = 0;
    hashdesc_init_basic(result->hashdesc, 16);
    return result;
}
ht64_t * ht64_alloc_size(size_t numel)
{
    ht64_t* result;
    result = malloc(sizeof(ht64_t));

    result->nalloc = numel;
    result->keys = (keyptr_t*)calloc(numel, sizeof(keyptr_t));
    result->hashdesc = hashdesc_alloc();
    result->nset = 0;
    hashdesc_init_basic(result->hashdesc, numel);
    return result;
}
void ht64_free(ht64_t* target) {
    size_t i;
    if (target->nalloc) {
        for (i = 0;i < target->nalloc;i++) {
            if (target->keys[i])
                _htkeydata_free64(target->keys[i]);
        }
        free(target->keys);
        target->nalloc = 0;
    }
    hashdesc_free(target->hashdesc);
    free(target);
}
int64_t ht64_get(ht64_t* target, char* key, size_t keylen, int* nullflag) {
    size_t hash_val;
    keyptr_t keyptr;
    int recorded_key_type;
    size_t recorded_key_length;
    int64_t result;
    int found;
    *nullflag = 1;
    if (!target)return 0;
    result = 0;
    hash_val = hash_index(key, keylen, target->hashdesc);
    if (target->keys[hash_val]) {
        keyptr = target->keys[hash_val];
        recorded_key_type = _htkeydata_type(keyptr);
        recorded_key_length = _htkeydata_length(keyptr);
        if (recorded_key_length == keylen && memcmp(key, _htkeydata_key(keyptr), keylen) == 0)
            found = 1;
        else
            found = 0;
        while (recorded_key_type == 2 && !found) {
            keyptr = _htkeydata_nextkey(keyptr, sizeof(int64_t));
            recorded_key_type = _htkeydata_type(keyptr);
            recorded_key_length = _htkeydata_length(keyptr);
            if (recorded_key_length == keylen && memcmp(key, _htkeydata_key(keyptr), keylen) == 0)
                found = 1;
        }
        if (found) {
            *nullflag = 0;
            result = _htkeydata_value64(keyptr);
        }
    }
    return result;
}
static void _ht64_rehash(ht64_t* target, size_t newnumel) {
    hashdesc_t* newhdesc;
    size_t i;
    keyptr_t* newkeys;
    keyptr_t tmpkey;
    keyptr_t tmpkey2;
    size_t keylen;
    int64_t keyvalue;
    size_t newindex;
    char* keyname;
    newhdesc = target->hashdesc;
    hashdesc_init_basic(newhdesc, newnumel);
    newkeys = (keyptr_t*)calloc(newnumel, sizeof(keyptr_t));
    for (i = 0;i < target->nalloc;i++) {
        if (target->keys[i]) {
            tmpkey = target->keys[i];
            while (_htkeydata_type(tmpkey) == 2)
            {
                keylen = _htkeydata_length(tmpkey);
                keyname = _htkeydata_key(tmpkey);
                keyvalue = _htkeydata_value64(tmpkey);
                newindex = hash_index((void*)keyname, keylen, newhdesc);
                if (newkeys[newindex] == NULL) {
                    newkeys[newindex] = tmpkey;
                    ((char*)tmpkey)[0] = 1; /* change the type to be an end node*/
                    tmpkey = _htkeydata_nextkey(tmpkey, sizeof(int64_t));
                }
                else {
                    tmpkey2 = newkeys[newindex];
                    while (_htkeydata_type(tmpkey2) == 2) {
                        tmpkey2 = _htkeydata_nextkey(tmpkey2,sizeof(int64_t));
                    }
                    _htkeydata_setnext(tmpkey2, tmpkey, sizeof(int64_t));
                    tmpkey = _htkeydata_nextkey(tmpkey, sizeof(int64_t));
                }
            }
            keylen = _htkeydata_length(tmpkey);
            keyname = _htkeydata_key(tmpkey);
            keyvalue = _htkeydata_value64(tmpkey);
            newindex = hash_index((void*)keyname, keylen, newhdesc);
            if (newkeys[newindex] == NULL) {
                newkeys[newindex] = tmpkey;
                ((char*)tmpkey)[0] = 1; /* change the type to be an end node*/
            }
            else {
                tmpkey2 = newkeys[newindex];
                while (_htkeydata_type(tmpkey2) == 2) {
                    tmpkey2 = _htkeydata_nextkey(tmpkey2, sizeof(int64_t));
                }
                _htkeydata_setnext(tmpkey2, tmpkey, sizeof(int64_t));
                tmpkey = _htkeydata_nextkey(tmpkey, sizeof(int64_t));
            }
        }
    }
    free(target->keys);
    target->keys = newkeys;
    target->nalloc = newnumel;
}
void ht64_set(ht64_t* target, char* key, size_t keylen, int64_t value, int *nullflag) {
    size_t hash_val;
    keyptr_t keyptr;
    size_t depth;
    *nullflag = 1;
    if (!target)return;
    hash_val = hash_index(key, keylen, target->hashdesc);
    if (target->keys[hash_val]) {
        keyptr = target->keys[hash_val];
        target->keys[hash_val] = _htkeydata_set64(keyptr, key, keylen, &depth, value, nullflag);
        if (*nullflag) {
            target->nset++;
            if (target->nset > target->nalloc)
                _ht64_rehash(target, target->nalloc * 2);
        }
    }
    else {
        target->keys[hash_val] = _htkeydata_create64(key, keylen, value);
        *nullflag = 1;
        target->nset++;
    }
}
void ht64_inc(ht64_t * target, char * key, size_t keylen, int64_t by, int * nullflag)
{
    size_t hash_val;
    keyptr_t keyptr;
    size_t depth;
    int64_t oldval;
    *nullflag = 1;
    if (!target)return;
    hash_val = hash_index(key, keylen, target->hashdesc);
    if (target->keys[hash_val]) {
        keyptr = target->keys[hash_val];
        oldval = _htkeydata_value64(keyptr);
        target->keys[hash_val] = _htkeydata_set64(keyptr, key, keylen, &depth, oldval + by, nullflag);
        if (*nullflag) {
            target->nset++;
            if (target->nset > target->nalloc)
                _ht64_rehash(target, target->nalloc * 2);
        }
    }
    else {
        target->keys[hash_val] = _htkeydata_create64(key, keylen, by);
        *nullflag = 1;
        target->nset++;
    }
}
size_t ht64_astables(ht64_t * target, char*** p_listofnames, size_t** p_listoflens, int64_t** p_values)
{
    char** result;
    size_t* resultlen;
    int64_t* resultvalues;
    size_t i, j;
    size_t len;
    keyptr_t curkey;

    if (p_values)
        resultvalues = (int64_t*)malloc(target->nset * sizeof(int64_t));
    else
        resultvalues = NULL;
    if (p_listofnames)
        result = (char**)malloc(target->nset * sizeof(char*));
    else
        result = NULL;
    if (p_listoflens)
        resultlen = (size_t*)malloc(target->nset * sizeof(size_t));
    else
        resultlen = NULL;
    j = 0;
    curkey = target->keys[0];
    for (i = 0;i < target->nset;i++) {
        if (target->keys[j] == NULL) {
            while (target->keys[j] == NULL)
                j++;
            curkey = target->keys[j];
        }
        len = _htkeydata_length(curkey);
        if (result) {
            result[i] = malloc((len + 1) * sizeof(char));
            memcpy(result[i], _htkeydata_key(curkey), len);
            result[i][len] = 0;
        }
        if (resultlen)
            resultlen[i] = len;
        if (resultvalues)
            resultvalues[i] = _htkeydata_value64(curkey);
        if (_htkeydata_type(curkey) == 2)
            curkey = _htkeydata_nextkey(curkey, sizeof(int64_t));
        else {
            j++;
            curkey = target->keys[j];
        }
    }
    if (p_listofnames)
        *p_listofnames = result;
    if (p_values)
        *p_values = resultvalues;
    if (p_listoflens)
        *p_listoflens = resultlen;
    return i;
}
void ht64_save(ht64_t* table, PF_t* file) {
    size_t i, j;
    size_t len;
    keyptr_t curkey;
    j = 0;
    curkey = table->keys[0];
    PFputint64(file, (uint64_t)(table->nalloc));
    PFputint64(file, (uint64_t)(table->nset));
    hashdesc_save(table->hashdesc, file);
    for (i = 0;i < table->nset;i++) {
        if (table->keys[j] == NULL) {
            while (table->keys[j] == NULL)
                j++;
            curkey = table->keys[j];
        }
        PFputint64(file, (int64_t)j);
        len = _htkeydata_length(curkey);
        PFputint64(file, (int64_t)len);
        PFwrite(_htkeydata_key(curkey), 1, len, file);
        PFputint64(file, _htkeydata_value64(curkey));
        if (_htkeydata_type(curkey) == 2)
            curkey = _htkeydata_nextkey(curkey, sizeof(int64_t));
        else {
            j++;
            curkey = table->keys[j];
        }
    }
}
ht64_t* ht64_load(PF_t* file) {
    ht64_t* result;
    size_t i,j;
    size_t len;
    char* key;
    keyptr_t keyptr;
    size_t depth;
    int nf;
    int64_t value;
    result = malloc(sizeof(ht64_t));
    result->nalloc = (size_t)PFgetint64(file);
    result->nset = (size_t)PFgetint64(file);
    result->hashdesc = hashdesc_load(file);
    result->keys = (keyptr_t*)calloc(result->nalloc, sizeof(keyptr_t));
    key = NULL;
    for (i = 0;i < result->nset;i++) {
        j = (size_t)PFgetint64(file);
        len = (size_t)PFgetint64(file);
        key = realloc(key,len);
        PFread(key, 1, len, file);
        value = PFgetint64(file);
        if (result->keys[j]) {
            keyptr = result->keys[j];
            result->keys[j] = _htkeydata_set64(keyptr, key, len, &depth, value, &nf);
        }
        else {
            result->keys[j] = _htkeydata_create64(key, len, value);
        }
    }
    if (key) free(key);
    return result;
}
