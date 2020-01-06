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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "textparsing.h"

/* useful string operations that bypass array preallocation */
char* text_append(char* target, char* to_append, size_t* len_after_start) {
    size_t delta;
#ifdef _WIN32
    strcat_s(target, *len_after_start, to_append);
#else
    strcat(target, to_append);
#endif
    delta = strlen(target);
    *len_after_start = *len_after_start - delta;
    return target + delta;
}
char* text_join(char** strings, char* separator, char* prefix, char* suffix, size_t count) {
    char* result;
    char* curstrptr;
    char* eosptr;
    size_t i;
    size_t totlen;
    size_t preflen;
    size_t seplen;
    size_t suflen;
    size_t remlen;

    preflen = 0;
    suflen = 0;
    seplen = 0;

    if (prefix) preflen = strlen(prefix);
    else prefix = "";
    if (separator) seplen = strlen(separator);
    else separator = "";
    if (suffix) suflen = strlen(suffix);
    else suffix = "";

    if (count == 0) {
        remlen = preflen + suflen + 1;
        result = (char*)calloc(preflen + suflen + 1, 1);
        if (preflen > 0) text_append(result, prefix, &remlen);
        if (suflen > 0) text_append(result, suffix, &remlen);
        return result;
    }

    totlen = preflen;
    for (i = 0;i < count;i++) {
        if (strings[i]) {
            totlen += strlen(strings[i]);
        }
    }
    totlen += suflen;
    totlen += seplen*(count - 1);
    result = (char*)calloc(totlen + 1, 1);
    remlen = totlen + 1;

    if (preflen > 0) eosptr = text_append(result, prefix, &remlen);
    else eosptr = result;
    i = 0;
    while (i < count - 1) {
        curstrptr = strings[i];
        eosptr = text_append(eosptr, curstrptr, &remlen);
        eosptr = text_append(eosptr, separator, &remlen);
        i++;
    }
    curstrptr = strings[i];
    eosptr = text_append(eosptr, curstrptr, &remlen);
    if (suflen > 0) text_append(eosptr, suffix, &remlen);
    return result;
}
char** text_splitshort(char* string, char* separator, size_t* outcount) {
    size_t i, j, k, delta;
    size_t nseps;
    char** result;
    char* start;
    i = 0;
    nseps = 0;
    j = 0;
    while (string[i]) {
        if (string[i] == separator[j]) j++;
        else j = 0;
        i++;
        if (!separator[j]) {
            nseps++;
            j = 0;
        }
    }
    result = malloc(sizeof(char*)*(nseps + 1));
    i = 0;
    j = 0;
    delta = 0;
    k = 0;
    start = string;
    while (string[i]) {
        if (string[i] == separator[j]) j++;
        else {
            delta += 1 + j;
            j = 0;
        }
        i++;
        if (!separator[j]) {
            result[k] = malloc(delta + 1);
            memcpy(result[k], start, delta);
            result[k][delta] = 0;
            start = string + i;
            k++;
            j = 0;
            delta = 0;
        }
    }
    result[k] = malloc(delta + 1);
    memcpy(result[k], start, delta);
    result[k][delta] = 0;
    *outcount = nseps + 1;
    return result;
}
char* text_fromint(int64_t value) {
    size_t length;
    char* result;
    if (value == 0) {
        result = calloc(2, sizeof(char));
        result[0] = '0';
        return result;
    }
    if (value > 0) {
        length = (size_t)log10((double)value) + 2;
    }
    else {
        length = (size_t)log10(0.0 - (double)value) + 3;
    }
    result = calloc(length, sizeof(char));
#ifdef _WIN32
    sprintf_s(result, length, _LLD_, (int64_t)value);
#else
    sprintf(result, _LLD_, (int64_t)value);
#endif
    return result;
}
char* text_fromdbl(double value, size_t digits) {
    size_t length;
    size_t extralength;
    size_t digits0;
    char* fmtspec;
    char* result;
    double logvalue;
    if (value == 0 || digits == 0) {
        result = calloc(2, sizeof(char));
        result[0] = '0';
        return result;
    }
    digits0 = digits;
    digits--;
    extralength = (size_t)(log10((double)(digits))) + 1;
    if (value > 0.0) logvalue = log10(value);
    else logvalue = log10(-value);
    if (value > 0.0) length = digits + 4;
    else length = digits + 5;

    fmtspec = calloc(extralength + 5, sizeof(char));
    fmtspec[0] = '%';
    fmtspec[1] = '.';
    if (logvalue >= (double)digits0) {
        length += (size_t)log10(logvalue) + 3;
    }
    else if (logvalue < 1.0 - (double)digits0) {
        length += (size_t)log10(-logvalue) + 3;
    }
    else {
        if (logvalue >= 1.0) digits -= (size_t)logvalue;
    }
#ifdef _WIN32
    sprintf_s(fmtspec + 2, extralength + 1, _LLD_, (int64_t)digits - 1);
#else
    sprintf(fmtspec + 2, _LLD_, (int64_t)digits - 1);
#endif
    if (logvalue >= (double)digits0) fmtspec[extralength + 2] = 'e';
    else if (logvalue < 1.0 - (double)digits0) fmtspec[extralength + 2] = 'e';
    else fmtspec[extralength + 2] = 'f';

    result = calloc(length, sizeof(char));
#ifdef _WIN32
    sprintf_s(result, length, fmtspec, value);
#else
    sprintf(result, fmtspec, value);
#endif
    free(fmtspec);
    return result;
}
size_t text_findchar(char* string, char c) {
    size_t i;
    i = 0;
    while (string[i] && string[i] != c)i++;
    return i;
}
size_t text_findstr(char* text, char* query) {
    char* result;
    result = strstr(text, query);
    if (result) return (size_t)(result - text);
    else return strlen(text);
}
size_t text_nextchar(const char* str, size_t start, char c) {
    while (str[start] && str[start]!=c)start++;
    return start;
}
size_t text_nextmatch(const char* str, size_t start, charcompare test) {
    while (str[start] && !(test(str[start])))start++;
    return start;
}
size_t text_nextspace(const char* str, size_t start) {
    while (str[start] && !((str[start] >= '\t' && str[start] <= '\r') || str[start] == ' ')) start++;
    return start;
}
size_t text_nextupper(const char* str, size_t start) {
    while (str[start] && (str[start]<'A' || str[start]>'Z'))start++;
    return start;
}
int iswordsep(char c) {
    return (c <= '0' || (c >= ':' && c <= '?') || (c >= '[' && c <= '^') || c == '`' || c >= '{');
}
size_t text_nextwordsep(const char* str, size_t start) {
    while (str[start] && !iswordsep(str[start]))start++;
    return start;
}
size_t text_nextdigit(const char* str, size_t start) {
    while (str[start] && str[start] >= '0' && str[start] <= '9') start++;
    return start;
}
int endswith(const char* end, const char* str) {
    size_t i;
    size_t j;
    i = strlen(str);
    j = strlen(end);
    while (i > 0 && j > 0) {
        --i; --j;
        if (str[i] != end[j])return 0;
    }
    if (i == 0 && j > 0)return 0;
    return 1;
}
int beginswith(const char* begin, const char* str) {
    size_t i, maxi;
    size_t j, maxj;
    maxi = strlen(str);
    maxj = strlen(begin);
    if (maxj > maxi)return 0;
    i = j = 0;
    while (i < maxi && j < maxj) {
        if (str[i] != begin[j])return 0;
        ++i; ++j;
    }
    if (i == maxi && j > maxj)return 0;
    return 1;
}