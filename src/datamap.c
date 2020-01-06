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

#include "datamap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct _descriptor_t{
    int algorithm;
    int unsorted;
} _descriptor_t;
struct _datamap32{
    int32_t* values;
    char** keys;
    int* keylens;
    size_t nelts;
    _descriptor_t descriptor;
    size_t allocated;
};
struct _datamap64{
    int64_t* values;
    char** keys;
    int* keylens;
    size_t nelts;
    _descriptor_t descriptor;
    size_t allocated;
};

DM32_t* new_DM32(int algorithm, int32_t default_value){
    DM32_t* result = NULL;
    if(algorithm==DM_ALGORITHM_BASICSORTEDLIST){
        result = (DM32_t*)malloc(sizeof(DM32_t));
        result->keys = NULL;
        result->values = NULL;
        result->keylens = NULL;
        result->descriptor.algorithm = DM_ALGORITHM_BASICSORTEDLIST;
        result->descriptor.unsorted = 1;
        result->nelts=0;
        result->allocated = 0;
    }
#ifdef CONSOLE_MODE
    else{
        fprintf(stderr,"The requested algorithm has not been implemented\n");
    }
#endif
    return result;
}
DM64_t* new_DM64(int algorithm, int64_t default_value){
    DM64_t* result = NULL;
    if(algorithm==DM_ALGORITHM_BASICSORTEDLIST){
        result = (DM64_t*)malloc(sizeof(DM64_t));
        result->keys = NULL;
        result->values = NULL;
        result->keylens = NULL;
        result->descriptor.algorithm = DM_ALGORITHM_BASICSORTEDLIST;
        result->descriptor.unsorted = 1;
        result->nelts=0;
        result->allocated = 0;
    }
    #ifdef CONSOLE_MODE
    else{
        fprintf(stderr,"The requested algorithm has not been implemented\n");
    }
    #endif
    return result;
}

void DM32_clear(DM32_t* target) {
    size_t celt;
    if (target->nelts == 0) {
        if (target->values) free(target->values);
        target->values = NULL;
        if (target->keys) free(target->keys);
        target->keys = NULL;
    }
    else {
        for (celt = 0;celt<target->nelts;celt++) {
            if (target->keys[celt]) free(target->keys[celt]);
            target->keys[celt] = 0;
        }
        free(target->values);
        target->values = NULL;
        free(target->keys);
        target->keys = NULL;
        free(target->keylens);
        target->keylens = NULL;
    }
    target->values = NULL;
    target->keys = NULL;
    target->nelts = 0;
    target->keylens = NULL;
}
void DM64_clear(DM64_t* target) {
    size_t celt;
    if (target->nelts == 0) {
        if (target->values) free(target->values);
        target->values = NULL;
        if (target->keys) free(target->keys);
        target->keys = NULL;
    }
    else {
        for (celt = 0;celt<target->nelts;celt++) {
            if (target->keys[celt]) free(target->keys[celt]);
            target->keys[celt] = 0;
        }
        free(target->values);
        target->values = NULL;
        free(target->keys);
        target->keys = NULL;
        free(target->keylens);
        target->keylens = NULL;
    }
    target->values = NULL;
    target->keys = NULL;
    target->nelts = 0;
    target->keylens = NULL;
}

void free_DM32(DM32_t* target){
    DM32_clear(target);
    free(target);
}
void free_DM64(DM64_t* target){
    DM64_clear(target);
    free(target);
}

int DM32_write(DM32_t* target, FILE* f) {
    size_t i;
    int64_t i64;
    int32_t i32;
    if (!target)return 1;
    DM32_sort(target);
    i32 = 0x1000; /* this is version 0.1.0, encoded as follows: .. ... ... */
    fwrite(&i32, sizeof(int32_t), 1, f);
    i64 = (int64_t) target->nelts;
    fwrite(&i64, sizeof(int64_t), 1, f);
    i32 = (int32_t) target->descriptor.algorithm;
    fwrite(&i32, sizeof(int32_t), 1, f);
    fwrite(target->values, sizeof(int32_t), target->nelts, f);
    fwrite(target->keylens, sizeof(int), target->nelts, f);
    for (i = 0;i < target->nelts;i++) {
        fwrite(target->keys[i], sizeof(char), target->keylens[i], f);
    }
    return 0;
}
int DM64_write(DM64_t* target, FILE* f) {
    size_t i;
    int64_t i64;
    int32_t i32;
    if (!target)return 1;
    DM64_sort(target);
    i32 = 0x1000; /* this is version 0.1.0, encoded as follows: .. ... ... */
    fwrite(&i32, sizeof(int32_t), 1, f);
    i64 = (int64_t)target->nelts;
    fwrite(&i64, sizeof(int64_t), 1, f);
    i32 = (int32_t)target->descriptor.algorithm;
    fwrite(&i32, sizeof(int32_t), 1, f);
    fwrite(target->values, sizeof(int64_t), target->nelts, f);
    fwrite(target->keylens, sizeof(int), target->nelts, f);
    for (i = 0;i < target->nelts;i++) {
        fwrite(target->keys[i], sizeof(char), target->keylens[i], f);
    }
    return 0;
}

int DM32_read(DM32_t* target, FILE* f) {
    size_t i;
    int64_t i64;
    int32_t i32;
    if (!target)return 1;
    DM32_clear(target);
    fread(&i32, sizeof(int32_t), 1, f);
    /* this is version 0.1.0, encoded as follows: .. ... ... */
    fread(&i64, sizeof(int64_t), 1, f);
    target->nelts = (size_t)i64;
    target->allocated = target->nelts;
    fread(&i32, sizeof(int32_t), 1, f);
    target->descriptor.algorithm = (int)i32;
    target->values = malloc(sizeof(int32_t)* target->nelts);
    target->keylens = malloc(sizeof(int)* target->nelts);
    target->keys = malloc(sizeof(char*)* target->nelts);
    fread(target->values, sizeof(int32_t), target->nelts, f);
    fread(target->keylens, sizeof(int), target->nelts, f);
    for (i = 0;i < target->nelts;i++) {
        target->keys[i] = malloc(sizeof(char)*target->keylens[i]);
        fread(target->keys[i], sizeof(char), target->keylens[i], f);
    }
    target->descriptor.unsorted = 0;
    return 0;
}
int DM64_read(DM64_t* target, FILE* f) {
    size_t i;
    int64_t i64;
    int32_t i32;
    if (!target)return 1;
    DM64_clear(target);
    fread(&i32, sizeof(int32_t), 1, f);
    /* this is version 0.1.0, encoded as follows: .. ... ... */
    fread(&i64, sizeof(int64_t), 1, f);
    target->nelts = (size_t)i64;
    target->allocated = target->nelts;
    fread(&i32, sizeof(int32_t), 1, f);
    target->descriptor.algorithm = (int)i32;
    target->values = malloc(sizeof(int64_t)* target->nelts);
    target->keylens = malloc(sizeof(int)* target->nelts);
    target->keys = malloc(sizeof(char*)* target->nelts);
    fread(target->values, sizeof(int64_t), target->nelts, f);
    fread(target->keylens, sizeof(int), target->nelts, f);
    for (i = 0;i < target->nelts;i++) {
        target->keys[i] = malloc(sizeof(char)*target->keylens[i]);
        fread(target->keys[i], sizeof(char), target->keylens[i], f);
    }
    target->descriptor.unsorted = 0;
    return 0;
}

int _tracking_cmpstr(const char* A, const char* B, size_t* begin, size_t Alen, size_t Blen){
    /* compare two string starting at (*begin) */
    while((*begin)<Alen && (*begin)<Blen && A[*begin]==B[*begin]){
        (*begin)++;
    }
    if((*begin)==Alen){
        if(Alen==Blen)
            return 0; /* A == B */
        else
            return -1; /* A < B */
    }
    if((*begin)==Blen)
        /* A > B  -- we know that A != B and (*begin) < min(Alen,Blen)*/
        return 1;
    if(A[*begin] < B[*begin])
        return -1;
    else
        return 1; /* A > B*/
}

size_t _DM32_findindexof_BASICSORTEDLIST(DM32_t* datamap, const char* key, int keylen, int* NULLflag){
    size_t minki, maxki, midki;
    size_t minbgn, maxbgn, midbgn;
    int maxcmp, mincmp, midcmp;

    *NULLflag = 0;

    minki = 0;
    maxki = datamap->nelts;
    minbgn = 0;
    maxbgn = 0;
    midbgn = 0;
    mincmp = _tracking_cmpstr(key,datamap->keys[minki],&minbgn,(size_t)keylen,(size_t)datamap->keylens[minki]);
    if(mincmp<0) *NULLflag=1;
    if(mincmp<=0) return minki;

    maxcmp = _tracking_cmpstr(key,datamap->keys[maxki-1],&maxbgn,(size_t)keylen,(size_t)datamap->keylens[maxki-1]);
    if(maxcmp>0){
        *NULLflag=1;
        return maxki;
    }
    if(maxcmp==0) return maxki-1;

    while(maxki - minki > 1){
        midbgn = (minbgn<maxbgn)? minbgn : maxbgn;
        midki = minki + (maxki-minki)/2;
        midcmp = _tracking_cmpstr(key,datamap->keys[midki],&midbgn,(size_t)keylen,(size_t)datamap->keylens[midki]);
        if(midcmp >= 0){
            minki=midki;
            minbgn=midbgn;
        }
        if(midcmp <= 0){
            maxki=midki;
            maxbgn=midbgn;
        }
    }
    if(minki!=maxki)
        *NULLflag=1;
    return maxki;
}
size_t _DM64_findindexof_BASICSORTEDLIST(DM64_t* datamap, const char* key, int keylen, int* NULLflag){
    size_t minki, maxki, midki;
    size_t minbgn, maxbgn, midbgn;
    int maxcmp, mincmp, midcmp;

    *NULLflag = 0;

    minki = 0;
    maxki = datamap->nelts;
    minbgn = 0;
    maxbgn = 0;
    midbgn = 0;
    mincmp = _tracking_cmpstr(key,datamap->keys[minki],&minbgn,(size_t)keylen,(size_t)datamap->keylens[minki]);
    if(mincmp<0) *NULLflag=1;
    if(mincmp<=0) return minki;

    maxcmp = _tracking_cmpstr(key,datamap->keys[maxki-1],&maxbgn,(size_t)keylen,(size_t)datamap->keylens[maxki-1]);
    if(maxcmp>0){
        *NULLflag=1;
        return maxki;
    }
    if(maxcmp==0) return maxki-1;

    while(maxki - minki > 1){
        midbgn = (minbgn<maxbgn)? minbgn : maxbgn;
        midki = minki + (maxki-minki)/2;
        midcmp = _tracking_cmpstr(key,datamap->keys[midki],&midbgn,(size_t)keylen,(size_t)datamap->keylens[midki]);
        if(midcmp >= 0){
            minki=midki;
            minbgn=midbgn;
        }
        if(midcmp <= 0){
            maxki=midki;
            maxbgn=midbgn;
        }
    }
    if(minki!=maxki)
        *NULLflag=1;
    return maxki;
}

void DM32_append(DM32_t* datamap, const char* key, int keylen, int32_t value){
    size_t indexv;
    if(datamap){
        datamap->descriptor.unsorted = 1;
        if(datamap->nelts == 0){
            if(datamap->values) free(datamap->values);
            if(datamap->keys) free(datamap->keys);
            /* we assume here that individual keys have properly been freed previously */
            datamap->values = (int32_t*)malloc(sizeof(int32_t));
            datamap->keys = (char**)malloc(sizeof(char*));
            datamap->keys[0] = (char*)malloc(sizeof(char)*keylen);
            datamap->keylens = (int*)malloc(sizeof(int));
            memcpy(datamap->keys[0],key,keylen);
            datamap->values[0] = value;
            datamap->nelts = 1;
            datamap->keylens[0] = keylen;
            datamap->allocated=1;
        }
        else{
            /* Create a new element regardless of whether one existed previously */
            if(datamap->descriptor.algorithm == DM_ALGORITHM_BASICSORTEDLIST){
                indexv = datamap->nelts;
                datamap->nelts++;
                if (datamap->nelts >= datamap->allocated) {
                    if (datamap->allocated*2 - datamap->nelts <= 0x1000000)
                        datamap->allocated *= 2;
                    else
                        datamap->allocated += 0x1000000;
                    datamap->values = (int32_t*)realloc(datamap->values, sizeof(int32_t)*datamap->allocated);
                    datamap->keylens = (int*)realloc(datamap->keylens, sizeof(int)*datamap->allocated);
                    datamap->keys = (char**)realloc(datamap->keys, sizeof(char*)*datamap->allocated);
                }
                datamap->values[indexv]=value;
                datamap->keylens[indexv]=keylen;
                datamap->keys[indexv] = (char*)malloc(sizeof(char)*keylen);
                memcpy(datamap->keys[indexv],key,keylen);
            }
        }
    }
}
void DM64_append(DM64_t* datamap, const char* key, int keylen, int64_t value){
    size_t indexv;
    if(datamap){
        datamap->descriptor.unsorted = 1;
        if(datamap->nelts == 0){
            if(datamap->values) free(datamap->values);
            if(datamap->keys) free(datamap->keys);
            /* we assume here that individual keys have properly been freed previously */
            datamap->values = (int64_t*)malloc(sizeof(int64_t));
            datamap->keys = (char**)malloc(sizeof(char*));
            datamap->keys[0] = (char*)malloc(sizeof(char)*keylen);
            datamap->keylens = (int*)malloc(sizeof(int));
            memcpy(datamap->keys[0],key,keylen);
            datamap->values[0] = value;
            datamap->nelts = 1;
            datamap->keylens[0] = keylen;
            datamap->allocated = 1;
        }
        else{
            /* Create a new element regardless of whether one existed previously */
            if (datamap->descriptor.algorithm == DM_ALGORITHM_BASICSORTEDLIST) {
                indexv = datamap->nelts;
                datamap->nelts++;
                if (datamap->nelts >= datamap->allocated) {
                    if(datamap->allocated*2 - datamap->nelts <= 0x100000000)
                        datamap->allocated *= 2;
                    else
                        datamap->allocated += 0x100000000;
                    datamap->values = (int64_t*)realloc(datamap->values, sizeof(int64_t)*datamap->allocated);
                    datamap->keylens = (int*)realloc(datamap->keylens, sizeof(int)*datamap->allocated);
                    datamap->keys = (char**)realloc(datamap->keys, sizeof(char*)*datamap->allocated);
                }
                datamap->values[indexv]=value;
                datamap->keylens[indexv]=keylen;
                datamap->keys[indexv] = (char*)malloc(sizeof(char)*keylen);
                memcpy(datamap->keys[indexv],key,keylen);
            }
        }
    }
}

void DM32_sort(DM32_t* datamap){
    size_t i1,i2,i3,i_n,i_n1,i_n2,starti,v1,v2;
    size_t tmpt;
    size_t *resultpartA;
    size_t *resultpartB;
    size_t *resultpart;
    size_t sort_level, new_sort_level;
    size_t newelts;
    int rel;
    int abortthis;

    DM32_t* tmpd;

    if(datamap && datamap->nelts>1 && datamap->descriptor.unsorted==1){
        if(datamap->descriptor.algorithm == DM_ALGORITHM_BASICSORTEDLIST){
            i1=0;
            i2=1;
            tmpt=0;
            /* based on mergesort.
            Since we expect that the list may be partially sorted, first run to the first
            unsorted element and begin from there. */
            while(i2<datamap->nelts && _tracking_cmpstr(datamap->keys[i1],datamap->keys[i2],&tmpt,datamap->keylens[i1],datamap->keylens[i2]) <= 0){
                tmpt = 0;
                i1++;
                i2++;
            }
            tmpt = 0;
            if(i2!=datamap->nelts){
                starti = i2;
                newelts = datamap->nelts-i2;
                resultpartA = (size_t*)malloc(sizeof(size_t)*newelts);
                resultpartB = (size_t*)malloc(sizeof(size_t)*newelts);
                resultpart = resultpartA;
                for(i_n=0;i_n<newelts;i_n++)
                    resultpart[i_n]=starti+i_n;
                sort_level = 1;
                while(sort_level<newelts){
                    /* outer loop: increasing sort level (log n)*/
                    i1=0;
                    i2=i1+sort_level;
                    new_sort_level = sort_level<<1;
                    while(i1<newelts){
                        /* inner loop 1: order all parts for a given sorting level (n/sort_level)*/
                        i_n1 = i1;
                        i_n2 = i2;
                        i3 = i2+sort_level;
                        if(i3>newelts) i3=newelts;
                        i_n = i1;
                        while(i_n<i3){
                            /* inner loop 2: order one block */
                            abortthis=0;
                            if(i_n1<i2)
                                v1 = resultpartA[i_n1];
                            else{
                                resultpartB[i_n] = resultpartA[i_n2];
                                i_n2++;
                                abortthis=1;
                            }
                            if(!abortthis){
                                if(i_n2<i3)
                                    v2 = resultpartA[i_n2];
                                else{
                                    resultpartB[i_n] = v1;
                                    i_n1++;
                                    abortthis=1;
                                }
                            }
                            if(!abortthis)
                            {
                                if(_tracking_cmpstr(datamap->keys[v1],datamap->keys[v2],&tmpt,datamap->keylens[v1],datamap->keylens[v2]) <= 0){
                                    resultpartB[i_n] = v1;
                                    i_n1++;
                                }
                                else{
                                    resultpartB[i_n] = v2;
                                    i_n2++;
                                }
                                tmpt = 0;

                            }
                            i_n++;
                        }
                        i1=i3;
                        i2=i1+sort_level;
                    }
                    /* the code below is necessary to avoid memory issues*/
                    resultpart=resultpartA;
                    resultpartA=resultpartB;
                    resultpartB=resultpart;
                    sort_level=new_sort_level;
                }
                free(resultpartB);resultpartB=NULL;
                resultpart = resultpartA;
                /* we now have resultpart as the sorted version of all previously unsorted elements */
                /* we need to find where to insert the old elements */
                resultpartB=(size_t*)malloc(sizeof(size_t)*datamap->nelts);
                i_n1=0;
                i_n2=0;
                i_n=0;
                while(i_n1<starti || i_n2<newelts){
                    /* inner loop 2: order one block */
                    abortthis=0;
                    if(i_n1<starti)
                        v1 = i_n1;
                    else{
                        resultpartB[i_n] = resultpartA[i_n2];
                        i_n2++;
                        abortthis=1;
                    }
                    if(!abortthis){
                        if(i_n2<newelts)
                            v2 = resultpartA[i_n2];
                        else{
                            resultpartB[i_n] = v1;
                            i_n1++;
                            abortthis=1;
                        }
                    }
                    if(!abortthis)
                    {
                        rel = _tracking_cmpstr(datamap->keys[v1],datamap->keys[v2],&tmpt,datamap->keylens[v1],datamap->keylens[v2]);
                        tmpt = 0;
                        if(rel <= 0){
                            resultpartB[i_n] = v1;
                            i_n1++;
                        }
                        else{
                            resultpartB[i_n] = v2;
                            i_n2++;
                        }

                    }
                    i_n++;
                }
                free(resultpart);
                resultpart=resultpartB;
                /* we now have resultpart as the sorted version all entries*/
                tmpd = new_DM32(DM_ALGORITHM_BASICSORTEDLIST, 0);
                tmpd->values = malloc(sizeof(int32_t)*datamap->allocated);
                tmpd->keys = malloc(sizeof(char*)*datamap->allocated);
                tmpd->keylens = malloc(sizeof(int)*datamap->allocated);
                tmpd->nelts = datamap->nelts;
                tmpd->allocated = datamap->allocated;
                memset(tmpd->keys,0,sizeof(char*)*datamap->nelts); /* NULL */
                i_n2=0;
                for(i_n=0;i_n<datamap->nelts;i_n++){
                    v1 = resultpart[i_n];
                    if(i_n2>0)
                        rel = _tracking_cmpstr(datamap->keys[v1],tmpd->keys[i_n2-1],&tmpt,datamap->keylens[v1],tmpd->keylens[i_n2-1]);
                    else
                        rel=1;
                    tmpt = 0;
                    if(rel==0){
                        /* we only keep the most recent entry, so the previous one needs to be free'd */
                        i_n2--;
                    }
                    if(rel<0){
                        /*shouldn't happen*/
                        fprintf(stderr,"Sorting the datamap failed!\n");
                    }
                    tmpd->values[i_n2] = datamap->values[v1];
                    tmpd->keys[i_n2] = datamap->keys[v1];
                    tmpd->keylens[i_n2] = datamap->keylens[v1];
                    i_n2++;
                }

                free(resultpart);
                free(datamap->values);
                free(datamap->keys);
                free(datamap->keylens);
                datamap->values = tmpd->values;
                datamap->keys = tmpd->keys;
                datamap->keylens = tmpd->keylens;
                datamap->nelts = i_n2;
                free(tmpd);
            }
        }
    }
    if(datamap){
        datamap->descriptor.unsorted=0;
    }
}
void DM64_sort(DM64_t* datamap){
    size_t i1,i2,i3,i_n,i_n1,i_n2,starti,v1,v2;
    size_t tmpt;
    size_t *resultpartA;
    size_t *resultpartB;
    size_t *resultpart;
    size_t sort_level, new_sort_level;
    size_t newelts;
    int rel;
    int abortthis;

    DM64_t* tmpd;

    if(datamap && datamap->nelts>1 && datamap->descriptor.unsorted == 1){
        if(datamap->descriptor.algorithm == DM_ALGORITHM_BASICSORTEDLIST){
            i1=0;
            i2=1;
            tmpt=0;
            /* based on mergesort.
            Since we expect that the list may be partially sorted, first run to the first
            unsorted element and begin from there. */
            while(i2<datamap->nelts && _tracking_cmpstr(datamap->keys[i1],datamap->keys[i2],&tmpt,datamap->keylens[i1],datamap->keylens[i2]) <= 0){
                tmpt = 0;
                i1++;
                i2++;
            }
            tmpt = 0;
            if(i2!=datamap->nelts){
                starti = i2;
                newelts = datamap->nelts-i2;
                resultpartA = (size_t*)malloc(sizeof(size_t)*newelts);
                resultpartB = (size_t*)malloc(sizeof(size_t)*newelts);
                resultpart = resultpartA;
                for(i_n=0;i_n<newelts;i_n++)
                    resultpart[i_n]=starti+i_n;
                sort_level = 1;
                while(sort_level<newelts){
                    /* outer loop: increasing sort level (log n)*/
                    i1=0;
                    i2=i1+sort_level;
                    new_sort_level = sort_level<<1;
                    while(i1<newelts){
                        /* inner loop 1: order all parts for a given sorting level (n/sort_level)*/
                        i_n1 = i1;
                        i_n2 = i2;
                        i3 = i2+sort_level;
                        if(i3>newelts) i3=newelts;
                        i_n = i1;
                        while(i_n<i3){
                            /* inner loop 2: order one block */
                            abortthis=0;
                            if(i_n1<i2)
                                v1 = resultpartA[i_n1];
                            else{
                                resultpartB[i_n] = resultpartA[i_n2];
                                i_n2++;
                                abortthis=1;
                            }
                            if(!abortthis){
                                if(i_n2<i3)
                                    v2 = resultpartA[i_n2];
                                else{
                                    resultpartB[i_n] = v1;
                                    i_n1++;
                                    abortthis=1;
                                }
                            }
                            if(!abortthis)
                            {
                                if(_tracking_cmpstr(datamap->keys[v1],datamap->keys[v2],&tmpt,datamap->keylens[v1],datamap->keylens[v2]) <= 0){
                                    resultpartB[i_n] = v1;
                                    i_n1++;
                                }
                                else{
                                    resultpartB[i_n] = v2;
                                    i_n2++;
                                }
                                tmpt = 0;

                            }
                            i_n++;
                        }
                        i1=i3;
                        i2=i1+sort_level;
                    }
                    /* the code below is necessary to avoid memory issues*/
                    resultpart=resultpartA;
                    resultpartA=resultpartB;
                    resultpartB=resultpart;
                    sort_level=new_sort_level;
                }
                free(resultpartB);resultpartB=NULL;
                resultpart = resultpartA;
                /* we now have resultpart as the sorted version of all previously unsorted elements */
                /* we need to find where to insert the old elements */
                resultpartB=(size_t*)malloc(sizeof(size_t)*datamap->nelts);
                i_n1=0;
                i_n2=0;
                i_n=0;
                while(i_n1<starti || i_n2<newelts){
                    abortthis=0;
                    if(i_n1<starti)
                        v1 = i_n1;
                    else{
                        /*all previously sorted elements have been added*/
                        resultpartB[i_n] = resultpartA[i_n2];
                        i_n2++;
                        abortthis=1;
                    }
                    if(!abortthis){
                        if(i_n2<newelts)
                            v2 = resultpartA[i_n2];
                        else{
                            /*all newly sorted elements have been added*/
                            resultpartB[i_n] = v1;
                            i_n1++;
                            abortthis=1;
                        }
                    }
                    if(!abortthis)
                    {
                        rel = _tracking_cmpstr(datamap->keys[v1],datamap->keys[v2],&tmpt,datamap->keylens[v1],datamap->keylens[v2]);
                        tmpt = 0;
                        if(rel <= 0){
                            resultpartB[i_n] = v1;
                            i_n1++;
                        }
                        else{
                            resultpartB[i_n] = v2;
                            i_n2++;
                        }
                    }
                    i_n++;
                }
                free(resultpart);
                resultpart=resultpartB;
                /* we now have resultpart as the sorted version all entries, but it might contain duplicates*/
                tmpd = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
                tmpd->values = malloc(sizeof(int64_t)*datamap->allocated);
                tmpd->keys = malloc(sizeof(char*)*datamap->allocated);
                tmpd->keylens = malloc(sizeof(int)*datamap->allocated);
                tmpd->nelts = datamap->nelts;
                tmpd->allocated = datamap->allocated;
                memset(tmpd->keys,0,sizeof(char*)*datamap->nelts); /* NULL */
                i_n2=0;
                for(i_n=0;i_n<datamap->nelts;i_n++){
                    v1 = resultpart[i_n];
                    if(i_n2>0)
                        rel = _tracking_cmpstr(datamap->keys[v1],tmpd->keys[i_n2-1],&tmpt,datamap->keylens[v1],tmpd->keylens[i_n2-1]);
                    else
                        rel=1;
                    tmpt = 0;
                    if(rel==0){
                        /* we only keep the most recent entry, so the previous one needs to be free'd */
                        i_n2--;
                    }
                    if(rel<0){
                        /*shouldn't happen*/
                        fprintf(stderr,"Sorting the datamap failed!\n");
                    }
                    tmpd->values[i_n2] = datamap->values[v1];
                    tmpd->keys[i_n2] = datamap->keys[v1];
                    tmpd->keylens[i_n2] = datamap->keylens[v1];
                    i_n2++;
                }

                free(resultpart);
                free(datamap->values);
                free(datamap->keys);
                free(datamap->keylens);
                datamap->values = tmpd->values;
                datamap->keys = tmpd->keys;
                datamap->keylens = tmpd->keylens;
                datamap->nelts = i_n2;
                free(tmpd);
            }
        }
    }
    if(datamap){
        datamap->descriptor.unsorted=0;
    }
}

void DM32_assign(DM32_t* datamap, const char* key, int keylen, int32_t value){
    int32_t* target_slot = NULL;
    int nullslot;
    size_t indexv;
    size_t i;
    if(keylen<=0) keylen = (int)strlen(key);
    if(datamap && datamap->descriptor.unsorted == 1) DM32_sort(datamap);
    if(datamap){
        if(datamap->nelts == 0){
            if(datamap->values) free(datamap->values);
            if(datamap->keys) free(datamap->keys);
            /* we assume here that individual keys have properly been freed previously */
            datamap->values = (int32_t*)malloc(sizeof(int32_t));
            datamap->keys = (char**)malloc(sizeof(char*));
            datamap->keys[0] = (char*)malloc(sizeof(char)*keylen);
            datamap->keylens = (int*)malloc(sizeof(int));
            memcpy(datamap->keys[0],key,keylen);
            datamap->values[0] = value;
            datamap->nelts = 1;
            datamap->keylens[0] = keylen;
            datamap->allocated = 1;
        }
        else{
            /* try and find an existing element, and create a new one if this fails. */
            /* this is where different algorithms come in play */
            if(datamap->descriptor.algorithm == DM_ALGORITHM_BASICSORTEDLIST){
                indexv = _DM32_findindexof_BASICSORTEDLIST(datamap,key,keylen,&nullslot);
                if(nullslot){
                    /* not found - reallocate & push everything forwards */
                    datamap->nelts++;
                    if (datamap->nelts >= datamap->allocated) {
                        if (datamap->allocated*2 - datamap->nelts <= 0x1000000)
                            datamap->allocated *= 2;
                        else
                            datamap->allocated += 0x1000000;
                        datamap->values = (int32_t*)realloc(datamap->values, sizeof(int32_t)*datamap->allocated);
                        datamap->keylens = (int*)realloc(datamap->keylens, sizeof(int)*datamap->allocated);
                        datamap->keys = (char**)realloc(datamap->keys, sizeof(char*)*datamap->allocated);
                    }
                    for(i=datamap->nelts-1;i>indexv;i--){
                        datamap->values[i]=datamap->values[i-1];
                        datamap->keylens[i]=datamap->keylens[i-1];
                        datamap->keys[i]=datamap->keys[i-1];
                    }
                }
                /* replace entry */
                datamap->values[indexv]=value;
                datamap->keylens[indexv]=keylen;
                if(!nullslot){
                    free(datamap->keys[indexv]);
                    datamap->keys[indexv] = NULL;
                }
                datamap->keys[indexv] = (char*)malloc(sizeof(char)*keylen);
                memcpy(datamap->keys[indexv],key,keylen);
            }
        }
    }
}
void DM64_assign(DM64_t* datamap, const char* key, int keylen, int64_t value){
    int64_t* target_slot = NULL;
    int nullslot;
    size_t indexv;
    size_t i;
    if(keylen<=0) keylen = (int)strlen(key);
    if(datamap){
        if(datamap->nelts == 0){
            if(datamap->values) free(datamap->values);
            if(datamap->keys) free(datamap->keys);
            /* we assume here that individual keys have properly been freed previously */
            datamap->values = (int64_t*)malloc(sizeof(int64_t));
            datamap->keys = (char**)malloc(sizeof(char*));
            datamap->keys[0] = (char*)malloc(sizeof(char)*keylen);
            datamap->keylens = (int*)malloc(sizeof(int));
            memcpy(datamap->keys[0],key,keylen);
            datamap->values[0] = value;
            datamap->nelts = 1;
            datamap->keylens[0] = keylen;
            datamap->allocated = 1;
        }
        else{
            /* try and find an existing element, and create a new one if this fails. */
            /* this is where different algorithms come in play */
            if(datamap->descriptor.algorithm == DM_ALGORITHM_BASICSORTEDLIST){
                indexv = _DM64_findindexof_BASICSORTEDLIST(datamap,key,keylen,&nullslot);
                if(nullslot){
                    /* not found - reallocate & push everything forwards */
                    datamap->nelts++;
                    if (datamap->nelts >= datamap->allocated) {
                        if (datamap->allocated*2 - datamap->nelts <= 0x100000000)
                            datamap->allocated *= 2;
                        else
                            datamap->allocated += 0x100000000;
                        datamap->values = (int64_t*)realloc(datamap->values, sizeof(int64_t)*datamap->allocated);
                        datamap->keylens = (int*)realloc(datamap->keylens, sizeof(int)*datamap->allocated);
                        datamap->keys = (char**)realloc(datamap->keys, sizeof(char*)*datamap->allocated);
                    }
                    for(i=datamap->nelts-1;i>indexv;i--){
                        datamap->values[i]=datamap->values[i-1];
                        datamap->keylens[i]=datamap->keylens[i-1];
                        datamap->keys[i]=datamap->keys[i-1];
                    }
                }
                /* replace entry */
                datamap->values[indexv]=value;
                datamap->keylens[indexv]=keylen;
                if(!nullslot){
                    free(datamap->keys[indexv]);
                    datamap->keys[indexv] = NULL;
                }
                datamap->keys[indexv] = (char*)malloc(sizeof(char)*keylen);
                memcpy(datamap->keys[indexv],key,keylen);
            }
        }
    }
}

int32_t DM32_get(DM32_t* datamap, const char* key, int keylen, int* NULLflag){
    int32_t result = 0;
    size_t telt;
    *NULLflag = 1;
    if(datamap && datamap->descriptor.unsorted == 1) DM32_sort(datamap);
    if(datamap && datamap->nelts > 0){
        telt = _DM32_findindexof_BASICSORTEDLIST(datamap,key,keylen,NULLflag);
        if(*NULLflag==0)
            result = datamap->values[telt];
    }
    return result;
}
int64_t DM64_get(DM64_t* datamap, const char* key, int keylen, int* NULLflag){
    int64_t result = 0;
    size_t telt;
    *NULLflag = 1;
    if(datamap && datamap->descriptor.unsorted == 1) DM64_sort(datamap);
    if(datamap && datamap->nelts > 0){
        telt = _DM64_findindexof_BASICSORTEDLIST(datamap,key,keylen,NULLflag);
        if(*NULLflag==0)
            result = datamap->values[telt];
    }
    return result;
}

void DM32_printf(DM32_t* datamap){
    size_t i;
    int j;
    if(datamap){
        for(i=0;i<datamap->nelts;i++){
            putchar('(');
            putchar('"');
            for(j=0;j<datamap->keylens[i];j++){
                putchar(datamap->keys[i][j]);
            }
            putchar('"');
            putchar(':');
            printf("%d",datamap->values[i]);
            putchar(')');
            if(i<datamap->nelts-1)
                putchar(',');
        }
    }
    else{
        printf("(null)");
    }
}
void DM64_printf(DM64_t* datamap){
    size_t i;
    int j;
    if(datamap){
        for(i=0;i<datamap->nelts;i++){
            putchar('(');
            putchar('"');
            for(j=0;j<datamap->keylens[i];j++){
                putchar(datamap->keys[i][j]);
            }
            putchar('"');
            putchar(':');
            printf("%zd",datamap->values[i]);
            putchar(')');
            if(i<datamap->nelts-1)
                putchar(',');
        }
    }
    else{
        printf("(null)");
    }
}

size_t DM32_numel(DM32_t * datamap)
{
    if (!datamap)return 0;
    return datamap->nelts;
}
size_t DM64_numel(DM64_t * datamap)
{
    if (!datamap)return 0;
    return datamap->nelts;
}

void DM32_ordertovalues(DM32_t* datamap) {
    size_t i;
    for (i = 0;i < datamap->nelts;i++) {
        datamap->values[i] = (int32_t)i;
    }
}
void DM64_ordertovalues(DM64_t* datamap) {
    size_t i;
    for (i = 0;i < datamap->nelts;i++) {
        datamap->values[i] = (int64_t)i;
    }
}
