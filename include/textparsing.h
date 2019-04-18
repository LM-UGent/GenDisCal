#ifndef DEF_TEXTPARSING
#define DEF_TEXTPARSING

#include <stdint.h>
#include <string.h>

typedef int(*charcompare)(char);

size_t text_nextchar(const char* str, size_t start, char c);
size_t text_nextmatch(const char* str, size_t start, charcompare test);
size_t text_nextspace(const char* str, size_t start);
size_t text_nextupper(const char* str, size_t start);
size_t text_nextwordsep(const char* str, size_t start);
size_t text_nextdigit(const char* str, size_t start);
int endswith(const char* end, const char* str);
int beginswith(const char* begin, const char* str);

#endif

