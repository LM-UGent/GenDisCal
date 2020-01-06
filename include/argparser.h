#ifndef DEF_ARGPARSER
#define DEF_ARGPARSER

#include <stdio.h>
#include <stdarg.h>

typedef struct argparser_t args_t;

args_t* args_alloc();
void args_free(args_t* target);

int args_parse(args_t* target, int argc, char** argv);
/* templateargs should be defined as a comma-separated list of the following:
str     this argument will be kept as is
int     this argument will be converted to an integer using atoi
float   this argument will be converted to a double using atof

...     (use this pattern by itself) unlimited amount of unconverted arguments
*/
void args_add(args_t* target, const char* longname, char shortname, const char* templateargs);
void args_add_help(args_t* target, const char* longname, const char* fullname, const char* desc, const char* shortdesc);

int args_getint(args_t* target, const char* longname, int argindex, int defaultval);
double args_getdouble(args_t* target, const char* longname, int argindex, double defaultval);
char* args_getstr(args_t* target, const char* longname, int argindex, char* defaultval);
int args_getintfromstr(args_t* target, const char* longname, int argindex, int defaultval);
double args_getdoublefromstr(args_t* target, const char* longname, int argindex, double defaultval);
int args_countargs(args_t* target, const char* longname);
int args_countfreeargs(args_t* target);
int args_ispresent(args_t* target, const char* longname);

void args_fprintgeneralhelp(FILE* target, args_t* arg, int colsperline);
void args_fprintspecifichelp(FILE* target, args_t* arg, const char* longname, int colsperline);
int args_get_quietness(args_t* target);
void args_report_error(args_t* target, const char* msg, ...);
void args_report_warning(args_t* target, const char* msg, ...);
void args_report_progress(args_t* target, const char* msg, ...);
void args_report_info(args_t* target, const char* msg, ...);
int args_is_helpmode(args_t* target);
/* warnlevels: VERBOSE,QUIETINFO,QUIETWARN,QUIETERR,QUIETALL*/
#define SUPERVERBOSE -1
#define VERBOSE 0
#define QUIETINFO 1
#define QUIETPROGRESS 2
#define QUIETWARN 3
#define QUIETERR 4
#define QUIETALL 5

#endif

