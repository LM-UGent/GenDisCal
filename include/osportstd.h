#ifndef DEF_OSPORTSTD
#define DEF_OSPORTSTD

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN32
#ifdef _WIN64
typedef int64_t ssize_t;
#else
typedef int32_t ssize_t;
#endif
#include <windows.h>
#include <io.h>
#else
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/wait.h>
#endif

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

#ifndef NAN
#define NAN (0/0)
#endif

#define rand_u32_minstd(x) (((uint64_t)x* 48271u) % 0x7fffffff)

typedef union any_t {
    uint64_t u64;
    uint32_t u32[2];
    uint16_t u16[4];
    uint8_t u8[8];
    int64_t i64;
    int32_t i32[2];
    int16_t i16[4];
    int8_t i8[8];
    double d;
    float f[2];
    char* str;
    void* ptr;
} any_t;

typedef int errcode_t;

int os_fileexists(const char* path);
char* os_rmdirname(char* path);
char* os_extractdirname(char* path);
char* os_stdoutfromexec(char* program, int argc, char** argv, size_t* outlen, size_t bufsize);
void* memcpyalloc(void* data, size_t datasize);

#define kilocount(x)    (x*1000)
#define Megacount(x)    (x*1000000)
#define Gigacount(x)    (x*1000000000)
#define Teracount(x)    (x*1000000000000)
#define Petacount(x)    (x*1000000000000000)
#define Exacount(x)     (x*1000000000000000000)
#define Zettacount(x)   (x*1000000000000000000000)
#define Yottacount(x)   (x*1000000000000000000000000)

#ifndef _LLD_
#ifdef _WIN32
#define _LLD_   "%lld"
#else
#define _LLD_   "%zd"
#endif
#endif


/*
Portable file handler
*/

typedef struct portablefile_t PF_t;

typedef char* unimem;

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
int PFseemslikeASCIItext(PF_t* f, size_t bytes_to_test);
char* PFgetfullpathptr(PF_t* f);
char* PFgetbasenameptr(PF_t* f);
FILE* PFgetFILE(PF_t* f);
size_t PFnumlines(PF_t* f, int countempty);

int16_t PFgetint16(PF_t* f);
int32_t PFgetint32(PF_t* f);
int64_t PFgetint64(PF_t* f);
float PFgetfloat(PF_t* f);
double PFgetdouble(PF_t* f);
/* note: PFreadline removes \n and \r at the end of lines */
char* PFreadline(PF_t* f);
unimem PFreadunimem16(PF_t* f);
unimem PFreadunimem32(PF_t* f);
unimem PFreadunimem64(PF_t* f);

errcode_t PFputint16(PF_t* f, int16_t input);
errcode_t PFputint32(PF_t* f, int32_t input);
errcode_t PFputint64(PF_t* f, int64_t input);
errcode_t PFputfloat(PF_t* f, float input);
errcode_t PFputdouble(PF_t* f, double input);

#endif
