#ifndef DEF_OSPORTSTD
#define DEF_OSPORTSTD

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#ifdef _WIN32
#ifdef _WIN64
typedef int64_t ssize_t;
#else
typedef int32_t ssize_t;
#endif
#include <io.h>
#else
#include <unistd.h>
#endif

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

#ifndef NAN
#define NAN (0/0)
#endif


typedef int errcode_t;

int os_fileexists(const char* path);
char* os_rmdirname(char* path);

#ifdef _WIN32
#define _LLD_   "%lld"
#else
#define _LLD_   "%zd"
#endif

/*
Portable file handler
*/

typedef struct portablefile_t PF_t;

#ifndef PFBUFFER
#define PFBUFFER 0x1000000
#endif

#ifndef PFSTDIODISABLE
#define PFSTDIOENABLE   1
#else
#define PFSTDIOENABLE   0
#endif

/* adapted std calls*/
errcode_t PFopen(PF_t** f, const char* path, const char* mode);
errcode_t PFclose(PF_t* f);
size_t PFread(void* buffer, size_t element_size, size_t element_count, PF_t* f);
size_t PFwrite(const void* buffer, size_t element_size, size_t element_count, PF_t* f);
long long PFtell(PF_t* f);
errcode_t PFseek(PF_t* f, long long offset, int origin);
errcode_t PFrewind(PF_t* f);
int PFgetc(PF_t* f);
errcode_t PFputc(int c, PF_t* f);
errcode_t PFprintf(PF_t* f, const char* fmt, ...);

/* extra functions */
FILE* PFgetFILE(PF_t* f);
size_t PFnumlines(PF_t* f, int countempty);

int16_t PFgetint16(PF_t* f);
int32_t PFgetint32(PF_t* f);
int64_t PFgetint64(PF_t* f);
float PFgetfloat(PF_t* f);
double PFgetdouble(PF_t* f);
char* PFreadline(PF_t* f);

errcode_t PFputint16(PF_t* f, int16_t input);
errcode_t PFputint32(PF_t* f, int32_t input);
errcode_t PFputint64(PF_t* f, int64_t input);
errcode_t PFputfloat(PF_t* f, float input);
errcode_t PFputdouble(PF_t* f, double input);

#endif
