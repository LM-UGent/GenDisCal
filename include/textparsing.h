#ifndef DEF_TEXTPARSING
#define DEF_TEXTPARSING

#include <stdint.h>
#include <string.h>

#ifndef _LLD_
#ifdef _WIN32
#define _LLD_   "%lld"
#else
#define _LLD_   "%zd"
#endif
#endif

typedef int(*charcompare)(char);

char* text_join(char** strings, char* separator, char* prefix, char* suffix, size_t count);
char** text_splitshort(char* string, char* separator, size_t* outcount);
char* text_fromint(int64_t value);
char* text_fromdbl(double value, size_t digits);
size_t text_findchar(char* string, char c);
size_t text_findstr(char* text, char* query);
size_t text_nextchar(const char* str, size_t start, char c);
size_t text_nextmatch(const char* str, size_t start, charcompare test);
size_t text_nextspace(const char* str, size_t start);
size_t text_nextupper(const char* str, size_t start);
int iswordsep(char c);
size_t text_nextwordsep(const char* str, size_t start);
size_t text_nextdigit(const char* str, size_t start);
int endswith(const char* end, const char* str);
int beginswith(const char* begin, const char* str);

#endif

