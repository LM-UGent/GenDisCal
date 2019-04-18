#include "taxextractor.h"
#include "datamap.h"
#include "textparsing.h"

struct taxextractor_t {
    DM64_t* names;
    int sorted;
};

taxextractor_t* taxextractor_alloc() {
    taxextractor_t* result;
    result = (taxextractor_t*) calloc(1, sizeof(taxextractor_t));
    result->names = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
    DM64_append(result->names, "\0", 1, 0);
    return result;
}
void taxextractor_free(taxextractor_t* target) {
    if (target) free_DM64(target->names);
    free(target);
}
char** _taxextractor_PCOF_S_to_str_array(const char* string, size_t* count) {
    size_t start, end, tmplen;
    int capcount;
    int nalloc;
    char** result;
    start = text_nextupper(string, 0);
    end = start;
    capcount = 0;
    nalloc = 16;
    result = calloc(16,sizeof(char*));
    while (string[end] != 0 && capcount<3) {
        end = text_nextupper(string, start + 1);
        tmplen = end - start;
        result[capcount] = realloc(result[capcount], tmplen + 6);
        result[capcount][0] = '0';
        result[capcount][1] = '0';
        result[capcount][2] = '0' + (char)capcount;
        result[capcount][3] = ':';
        result[capcount][4] = ':';
        memcpy(result[capcount] + 5, string + start, tmplen);
        result[capcount][tmplen + 5] = 0;
        capcount++;
        start = end;
    }
    while (string[end] != 0) {
        if (capcount >= nalloc) {
            nalloc *= 2;
            result = realloc(result,nalloc*sizeof(char*));
        }
        end = text_nextchar(string, start + 1, '_');
        tmplen = end - start;
        result[capcount] = realloc(result[capcount], tmplen + 6);
        result[capcount][0] = '0' + (char)((capcount / 100) % 10);
        result[capcount][1] = '0' + (char)((capcount / 10) % 10);
        result[capcount][2] = '0' + (char)(capcount % 10);
        result[capcount][3] = ':';
        result[capcount][4] = ':';
        memcpy(result[capcount] + 5, string + start, tmplen);
        result[capcount][tmplen + 5] = 0;
        capcount++;
        start = end + 1;
    }
    *count = capcount;
    result = realloc(result, capcount * sizeof(char*));
    return result;
}
void taxextractor_add_PCOF_S(taxextractor_t* extractor, const char* string) {
    char** tmpstr;
    size_t i, numel;
    extractor->sorted = 0;
    tmpstr = _taxextractor_PCOF_S_to_str_array(string, &numel);
    for (i = 0;i < numel;i++) {
        DM64_append(extractor->names, tmpstr[i], (int) strlen(tmpstr[i]), 0);
        free(tmpstr[i]);
    }
    free(tmpstr);
}

int taxextractor_translate_PCOF_S(taxextractor_t* extractor, const char* string, int64_t* idarray, size_t idarraysize) {
    char** tmpstr;
    size_t i, numel;
    int nullflag;
    if (!extractor->sorted) {
        DM64_sort(extractor->names);
        DM64_ordertovalues(extractor->names);
        extractor->sorted = 1;
    }
    tmpstr = _taxextractor_PCOF_S_to_str_array(string, &numel);
    for (i = 0;i < numel;i++) {
        if (i < idarraysize) {
            idarray[i] = DM64_get(extractor->names, tmpstr[i], (int)strlen(tmpstr[i]), &nullflag);
            if (nullflag)
                idarray[i] = 0;
        }
        free(tmpstr[i]);
    }
    while (i < idarraysize) {
        idarray[i] = 0;
        i++;
    }
    free(tmpstr);
    return 0;
}

int taxextractor_add(taxextractor_t* extractor, const char* method, const char* string) {
    if (strcmp(method, "PCOF_S") == 0) {
        taxextractor_add_PCOF_S(extractor, string);
        return 0;
    }
    return 1;
}

int taxextractor_translate(taxextractor_t* extractor, const char* method, const char* string, int64_t* idarray, size_t idarraysize) {
    if (strcmp(method, "PCOF_S") == 0) {
        return taxextractor_translate_PCOF_S(extractor, string, idarray, idarraysize);
    }
    return 1;
}

void taxextractor_sort(taxextractor_t* extractor) {
    if (!extractor->sorted) {
        DM64_sort(extractor->names);
        DM64_ordertovalues(extractor->names);
        extractor->sorted = 1;
    }
}