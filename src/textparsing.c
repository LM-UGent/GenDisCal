#include "textparsing.h"

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
inline int iswordsep(char c) {
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