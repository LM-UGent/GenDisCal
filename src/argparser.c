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

#include "argparser.h"
#include "datamap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define datamap_t DM32_t
#define free_datamap free_DM32
#define alloc_datamap new_DM32
#define datamap_assign DM32_assign
#define datamap_append DM32_append
#define datamap_sort DM32_sort
#define datamap_get DM32_get
#define datamap_int_t int32_t

typedef struct argument_t {
    char* longname;
    char shortname;
    int maxargcount;
    int* iargs;
    double* fargs;
    char** sargs;
    int niargs;
    int nfargs;
    int nsargs;
    short* iargset;
    short* fargset;
    short* sargset;
    char* argpattern;
    char* fullname;
    char* description;
    char* shortdesc;
    int cargvalue;
    int cargvaluei;
    int cargvaluef;
    int cargvalues;
    int argpresent;
    int occurencecount;
} argument_t;

argument_t* argument_alloc() {
    argument_t* result;
    result = (argument_t*)malloc(sizeof(argument_t));
    result->longname = NULL;
    result->shortname = 0;
    result->maxargcount = 0;
    result->iargs = NULL;
    result->fargs = NULL;
    result->sargs = NULL;
    result->iargset = NULL;
    result->fargset = NULL;
    result->sargset = NULL;
    result->niargs = 0;
    result->nfargs = 0;
    result->nsargs = 0;
    result->argpattern = NULL;
    result->fullname = NULL;
    result->description = NULL;
    result->shortdesc = NULL;
    result->cargvalue = 0;
    result->cargvaluei = 0;
    result->cargvaluef = 0;
    result->cargvalues = 0;
    result->argpresent = 0;
    result->occurencecount = 0;
    return result;
}

void argument_free(argument_t* target) {
    int i;
    if (target->longname)free(target->longname);
    if (target->iargs)free(target->iargs);
    if (target->fargs)free(target->fargs);
    if (target->sargs){
        for (i = 0;i < target->nsargs;i++) {
            if (target->sargs[i])free(target->sargs[i]);
        }
        free(target->sargs);
    }
    if (target->iargset)free(target->iargset);
    if (target->fargset)free(target->fargset);
    if (target->argpattern) free(target->argpattern);
    if (target->fullname) free(target->fullname);
    if (target->description) free(target->description);
    if (target->shortdesc) free(target->shortdesc);
    free(target);
}

void argument_setname(argument_t* target, const char* longver, char shortver) {
    size_t lverlen;
    if (target) {
        /* cleanup previous values */
        if (target->longname)free(target->longname);
        /* set new values */
        lverlen = strlen(longver);
        target->longname = (char*)malloc(lverlen + 1);
        memcpy(target->longname, longver, lverlen + 1);
        target->shortname = shortver;
    }
}
void argument_setargpattern(argument_t* target, const char* argpattern) {
    size_t curindex;
    size_t ci;
    int i;
    if (target) {
        /* cleanup previous values */
        if (target->iargs)free(target->iargs);
        if (target->fargs)free(target->fargs);
        if (target->sargs) {
            for (i = 0;i < target->nsargs;i++) {
                if (target->sargs[i])free(target->sargs[i]);
            }
            free(target->sargs);
        }
        if (target->argpattern) free(target->argpattern);
        /* set new values */
        if (strcmp(argpattern, "...") != 0) {
            curindex = 0;
            while (argpattern[curindex] != 0) {
                if (argpattern[curindex] == 'i' &&
                    argpattern[curindex + 1] == 'n' &&
                    argpattern[curindex + 2] == 't') {
                    target->niargs++;
                    target->maxargcount++;
                    curindex += 3;
                }
                else if (argpattern[curindex] == 'f' &&
                    argpattern[curindex + 1] == 'l' &&
                    argpattern[curindex + 2] == 'o' &&
                    argpattern[curindex + 3] == 'a' &&
                    argpattern[curindex + 4] == 't') {
                    target->nfargs++;
                    target->maxargcount++;
                    curindex += 5;
                }
                else if (argpattern[curindex] == 'v' &&
                    argpattern[curindex + 1] == 'a' &&
                    argpattern[curindex + 2] == 'r') {
                    target->nsargs++;
                    target->maxargcount++;
                    curindex += 3;
                }
                else if (argpattern[curindex] == 's' &&
                    argpattern[curindex + 1] == 't' &&
                    argpattern[curindex + 2] == 'r') {
                    target->nsargs++;
                    target->maxargcount++;
                    curindex += 3;
                }
                while (argpattern[curindex] != ',' && argpattern[curindex] != 0)
                    curindex++;
                if (argpattern[curindex] != 0)
                    curindex++;
            }
            target->argpattern = (char*)malloc(target->maxargcount + 1);
            if (target->niargs > 0) {
                target->iargs = (int*)malloc(sizeof(int)*target->niargs);
                target->iargset = (short*)malloc(sizeof(short)*target->niargs);
                memset(target->iargs, 0, sizeof(int)*target->niargs);
                memset(target->iargset, 0, sizeof(short)*target->niargs);
            }
            if (target->nfargs > 0) {
                target->fargs = (double*)malloc(sizeof(double)*target->nfargs);
                target->fargset = (short*)malloc(sizeof(short)*target->nfargs);
                memset(target->fargs, 0, sizeof(double)*target->nfargs);
                memset(target->fargset, 0, sizeof(short)*target->nfargs);
            }
            if (target->nsargs > 0) {
                target->sargs = (char**)malloc(sizeof(char*)*target->nsargs);
                target->sargset = (short*)malloc(sizeof(short)*target->nsargs);
                memset(target->sargs, 0, sizeof(char*)*target->nsargs);
                memset(target->sargset, 0, sizeof(short)*target->nsargs);
            }            
            memset(target->argpattern, 0, target->maxargcount + 1);
            curindex = 0;
            ci = 0;
            while (argpattern[curindex] != 0) {
                if (argpattern[curindex] == 'i' &&
                    argpattern[curindex + 1] == 'n' &&
                    argpattern[curindex + 2] == 't') {
                    target->argpattern[ci] = 'i';
                    ci++;
                }
                else if (argpattern[curindex] == 'f' &&
                    argpattern[curindex + 1] == 'l' &&
                    argpattern[curindex + 2] == 'o' &&
                    argpattern[curindex + 3] == 'a' &&
                    argpattern[curindex + 4] == 't') {
                    target->argpattern[ci] = 'f';
                    ci++;
                }
                else if (argpattern[curindex] == 's' &&
                    argpattern[curindex + 1] == 't' &&
                    argpattern[curindex + 2] == 'r') {
                    target->argpattern[ci] = 's';
                    ci++;
                }
                else if (argpattern[curindex] == 'v' &&
                    argpattern[curindex + 1] == 'a' &&
                    argpattern[curindex + 2] == 'r') {
                    target->argpattern[ci] = 'v';
                    ci++;
                }
                while (argpattern[curindex] != ',' && argpattern[curindex] != 0)
                    curindex++;
                if (argpattern[curindex] != 0)
                    curindex++;
            }
            target->argpattern[ci] = 0;
        }
        else {
            target->maxargcount = -1;
            target->niargs = 0;
            target->nfargs = 0;
            target->nsargs = 0;
            target->argpattern = (char*)malloc(2);
            target->argpattern[0] = '.';
            target->argpattern[1] = 0;
        }
    }
}
void argument_sethelp(argument_t* target, const char* fullname, const char* description, const char* shortdesc) {
    size_t fnlen;
    size_t dlen;
    size_t sdlen;
    if (target) {
        /* cleanup previous values */
        if (target->fullname) free(target->fullname);
        if (target->description) free(target->description);
        if (target->shortdesc) free(target->shortdesc);
        /* set new values */
        if (fullname) {
            fnlen = strlen(fullname) + 1;
            target->fullname = malloc(fnlen);
            memcpy(target->fullname, fullname, fnlen);
        }
        else
            target->fullname = NULL;
        if (description) {
            dlen = strlen(description) + 1;
            target->description = malloc(dlen);
            memcpy(target->description, description, dlen);
        }
        else
            target->description = NULL;
        if (shortdesc) {
            sdlen = strlen(shortdesc) + 1;
            target->shortdesc = malloc(sdlen);
            memcpy(target->shortdesc, shortdesc, sdlen);
        }
        else
            target->shortdesc = NULL;
    }
}

int argument_argstart(argument_t* target, int is_argv) {
    if (target) {
        if (target->maxargcount >= 0) {
            target->cargvalue = 0;
            target->cargvalues = 0;
        }
        target->cargvaluei = 0;
        target->cargvaluef = 0;
        target->argpresent = is_argv;
        if (is_argv)
            target->occurencecount++;
        if (target->maxargcount != 0)return 1;
    }
    return 0;
}
int argument_nextargvalue(argument_t* target, const char* value) {
    size_t valuelen;
    int output;
    /*Output values:
        0: argvalue properly added, but it's the last one
        1: argvalue properly added
        2: argvalue not added
    */
    if (target && (target->maxargcount < 0 || target->cargvalue < target->maxargcount)) {
        if (target->maxargcount < 0) {
            valuelen = strlen(value) + 1;
            target->sargs = (char**)realloc(target->sargs, sizeof(char*)*(target->nsargs + 1));
            target->sargs[target->nsargs] = (char*)malloc(valuelen);
            memcpy(target->sargs[target->nsargs], value, valuelen);
            target->nsargs++;
            target->cargvalues++;
            output = 1;
        }
        else {
            switch (target->argpattern[target->cargvalue]) {
            case 'i':
                target->iargs[target->cargvaluei] = atoi(value);
                target->iargset[target->cargvaluei] = 1;
                target->cargvaluei++;
                break;
            case 'f':
                target->fargs[target->cargvaluef] = atof(value);
                target->fargset[target->cargvaluef] = 1;
                target->cargvaluef++;
                break;
            case 's':
            case 'v':
                valuelen = strlen(value) + 1;
                target->sargs[target->cargvalues] = (char*)malloc(valuelen);
                target->sargset[target->cargvalues] = 1;
                memcpy(target->sargs[target->cargvalues], value, valuelen);
                target->cargvalues++;
                break;
            }
            target->cargvalue++;
            if (target->cargvalue == target->maxargcount) output = 0;
            else output = 1;
        }
    }
    else
        output = 2;
    return output;
}

int argument_geti(argument_t* target, int iid, int defaultval) {
    if (target && iid < target->niargs && target->iargset[iid]) {
        return target->iargs[iid];
    }
    else return defaultval;
}
double argument_getf(argument_t* target, int fid, double defaultval) {
    if (target && fid < target->nfargs && target->fargset[fid]) {
        return target->fargs[fid];
    }
    else return defaultval;
}
char* argument_gets(argument_t* target, int sid, char* defaultval) {
    if (target->sargset) {
        if (target && sid < target->nsargs && target->sargset[sid]) {
            return target->sargs[sid];
        }
    }
    else {
        if (target && sid < target->nsargs) {
            return target->sargs[sid];
        }
    }
    return defaultval;
}

struct argparser_t {
    datamap_t* argid;
    argument_t** argvs;
    size_t nargts;
    int quiet_level;
    int helpmodeon;
};

args_t* args_alloc(){
    args_t* result;
    result = malloc(sizeof(args_t));

    result->argid = alloc_datamap(DM_ALGORITHM_BASICSORTEDLIST, 0);
    result->argvs = (argument_t**)malloc(sizeof(argument_t*)*3);
    result->argvs[0] = argument_alloc();
    argument_setargpattern(result->argvs[0], "...");
    result->nargts = 3;
    /* Common commands for all programs: --help/-h and --quiet/-q*/
    result->argvs[1] = argument_alloc();
    datamap_assign(result->argid, "help", 4, 1);
    datamap_assign(result->argid, "h", 1, 1);
    argument_setname(result->argvs[1], "help", 'h');
    argument_setargpattern(result->argvs[1], "int");
    argument_sethelp(result->argvs[1], "HELP",
        "Displays help for individual arguments or everything if used alone. The integer parameter is used to define the printing width. (default: 80)",
        "Print help for specified arguments, or this help if none are specified");
    
    result->argvs[2] = argument_alloc();
    datamap_assign(result->argid, "quiet", 5, 2);
    datamap_assign(result->argid, "q", 1, 2);
    argument_setname(result->argvs[2], "quiet", 'q');
    argument_setargpattern(result->argvs[2], "str");
    argument_sethelp(result->argvs[2], "QUIET MODE",
        "Will mask some of the output of the program, depending on the argument specified.\nOptions are as follows:\n"
        "\n  none     No output will be masked, and internal steps will be shown."
        "\n  info     No information output will be produced."
        "\n  progress Only errors and warnings will be shown."
        "\n  warnings Only errors leading to program failure will be shown."
        "\n  all      No output will be produced\n",
        "quiet mode ");
    result->quiet_level = 1;
    result->helpmodeon = 1;
    return result;
}
void args_free(args_t* target){
    size_t i;
    free_datamap(target->argid);
    for (i = 0;i < target->nargts;i++) {
        argument_free(target->argvs[i]);
        target->argvs[i] = NULL;
    }
    free(target->argvs);
    target->argvs = NULL;
    free(target);
}

void args_summary_long(args_t* target,const char* prefix) {
    size_t tmpi, j;
    for (tmpi = 0;tmpi < target->nargts;tmpi++) {
        if (target->argvs[tmpi]->argpresent) {
            if (target->argvs[tmpi]->longname) {
                fprintf(stderr, "%s  X  %s \t\t", prefix, target->argvs[tmpi]->longname);
                fprintf(stderr, "(int x %d, ", target->argvs[tmpi]->niargs);
                fprintf(stderr, "float x %d, ", target->argvs[tmpi]->nfargs);
                fprintf(stderr, "string x %d)", target->argvs[tmpi]->nsargs);
                if (target->argvs[tmpi]->niargs > 0) {
                    fprintf(stderr, "\n%s           int: ", prefix);
                    for (j = 0;j < target->argvs[tmpi]->niargs;j++) {
                        fprintf(stderr, "%d, ", target->argvs[tmpi]->iargs[j]);
                    }
                }
                if (target->argvs[tmpi]->nfargs > 0) {
                    fprintf(stderr, "\n%s         float: ", prefix);
                    for (j = 0;j < target->argvs[tmpi]->nfargs;j++) {
                        fprintf(stderr, "%f, ", target->argvs[tmpi]->fargs[j]);
                    }
                }
                if (target->argvs[tmpi]->nsargs > 0) {
                    fprintf(stderr, "\n%s        string: ",prefix);
                    for (j = 0;j < target->argvs[tmpi]->nsargs;j++) {
                        fprintf(stderr, "\"%s\", ", target->argvs[tmpi]->sargs[j]);
                    }
                }
                fprintf(stderr, "\n");
            }
            else {
                fprintf(stderr, "%s  Unnamed arguments \t", prefix);
                fprintf(stderr, " x %d", target->argvs[tmpi]->nsargs);
                if (target->argvs[tmpi]->nsargs > 0)
                    fprintf(stderr, "\n%s        ", prefix);
                for (j = 0;j < target->argvs[tmpi]->nsargs;j++) {
                    fprintf(stderr, "\"%s\", ", target->argvs[tmpi]->sargs[j]);
                }
                fprintf(stderr, "\n");
            }

        }
        else {
            if (target->argvs[tmpi]->longname) {
                fprintf(stderr, "%s     %s\n", prefix, target->argvs[tmpi]->longname);
            }
            else {
                fprintf(stderr, "%s  Unnamed arguments\n",prefix);
            }
        }
    }
}

int args_parse(args_t* target, int argc, char** argv){
    int carg;
    int nullflag;
    size_t tmpi;
    size_t argtid;
    size_t colsperline;
    int parse_output;
    int untyped_argument_init=0;
    int anyargpresent;
    char* sargtype;
    char* tmpstr;
    carg = 1;
    argtid = 0;
    nullflag = 0;
    argument_argstart(target->argvs[argtid], 1);
    while (carg < argc) {
        if (argv[carg][0] == '-') {
            nullflag = 0;
            if (argv[carg][1] == '-') {
                /* long version */
                sargtype = argv[carg] + 2;
                argtid = (size_t)datamap_get(target->argid, sargtype, (int)strlen(sargtype), &nullflag);
                if (nullflag) {
                    fprintf(stderr, "[WARN] Invalid argument \"%s\" ignored\n", argv[carg]);
                    argtid = 0;
                }
                else {
                    parse_output = argument_argstart(target->argvs[argtid], 1);
                    if (parse_output == 0)
                        argtid = 0;
                }
            }
            else {
                /* short version */
                tmpi = 1;
                while (argv[carg][tmpi]) {
                    sargtype = argv[carg] + tmpi;
                    argtid = (size_t)datamap_get(target->argid, sargtype, 1, &nullflag);
                    if (nullflag) {
                        fprintf(stderr, "[WARN] Invalid argument \"-%c\" ignored\n", argv[carg][tmpi]);
                        argtid = 0;
                    }
                    else {
                        parse_output = argument_argstart(target->argvs[argtid], 1);
                        if (parse_output == 0)
                            argtid = 0;
                    }
                    tmpi++;
                }
            }
        }
        else {
            parse_output = argument_nextargvalue(target->argvs[argtid], argv[carg]);
            if (parse_output == 0 || parse_output == 2)
                argtid = 0;
        }
        carg++;
    }
    /* quiet level */
    if (target->argvs[2]->argpresent) {
        tmpstr = argument_gets(target->argvs[2], 0, NULL);
        if(!tmpstr){
            target->quiet_level = 3; /*errors will still be shown*/
        }
        else {
            if (strcmp("debug", tmpstr) == 0) target->quiet_level = SUPERVERBOSE;
            if (strcmp("none", tmpstr) == 0) target->quiet_level = VERBOSE;
            if (strcmp("info", tmpstr) == 0) target->quiet_level = QUIETINFO;
            if (strcmp("progress", tmpstr) == 0) target->quiet_level = QUIETWARN;
            if (strcmp("warnings", tmpstr) == 0) target->quiet_level = QUIETERR;
            if (strcmp("all", tmpstr) == 0) target->quiet_level = 4;
        }
    }
    args_report_info(target, NULL);
    args_report_progress(target, NULL);
    args_report_warning(target, NULL);
    args_report_error(target, NULL);

    /* help is shown? is yes, output 1, otherwise output 0*/
    if (target->argvs[1]->argpresent) {
        colsperline = args_getint(target, "help", 0, 80);
        anyargpresent = 0;
        fprintf(stderr, "\n");
        if (target->argvs[1]->argpresent && target->argvs[1]->occurencecount>1) {
            args_fprintspecifichelp(stderr, target, target->argvs[1]->longname, (int)colsperline);
            anyargpresent = 1;
        }
        for (tmpi = 2;tmpi < target->nargts;tmpi++) {
            if (target->argvs[tmpi]->argpresent) {
                if (anyargpresent)
                    fprintf(stderr, "\n");
                args_fprintspecifichelp(stderr, target, target->argvs[tmpi]->longname, (int)colsperline);
                anyargpresent = 1;
            }
        }
        if (anyargpresent == 0) {
            args_fprintgeneralhelp(stderr, target, (int)colsperline);
        }
        return 1;
    }
    else {
        target->helpmodeon = 0;
    }
    if (target->quiet_level == SUPERVERBOSE) {
        args_report_info(NULL, "List of arguments:\n");
        args_summary_long(target,"[INFO]");
    }
    return 0;
}
void args_add(args_t* target, const char* longname, char shortname, const char* templateargs){
    if (target) {
        target->argvs = realloc(target->argvs, sizeof(argument_t*)*(target->nargts+1));
        target->argvs[target->nargts] = argument_alloc();
        argument_setargpattern(target->argvs[target->nargts], templateargs);
        argument_setname(target->argvs[target->nargts], longname, shortname);
        datamap_assign(target->argid, longname, (int)strlen(longname), (datamap_int_t)target->nargts);
        datamap_assign(target->argid, &shortname, 1, (datamap_int_t)target->nargts);
        target->nargts++;
    }
}
void args_add_help(args_t* target, const char* longname, const char* fullname, const char* desc, const char* shortdesc) {
    size_t index;
    int nullflag;
    if (target) {
        index = (size_t)datamap_get(target->argid, longname, (int)strlen(longname), &nullflag);
        argument_sethelp(target->argvs[index], fullname, desc, shortdesc);
    }
}

int args_getint(args_t* target, const char* longname, int argindex, int defaultval){
    size_t index;
    int nullflag=0;
    if (target) {
        if (!longname) {
            index = 0;
            argument_geti(target->argvs[index], argindex, defaultval);
        }
        else {
            index = (size_t)datamap_get(target->argid, longname, (int)strlen(longname), &nullflag);
            if (!nullflag && target->argvs[index]->argpresent)
                return argument_geti(target->argvs[index], argindex, defaultval);
        }
    }
    return defaultval;
}
double args_getdouble(args_t* target, const char* longname, int argindex, double defaultval){
    size_t index;
    int nullflag=0;
    if (target) {
        if (!longname) {
            index = 0;
            argument_getf(target->argvs[index], argindex, defaultval);
        }
        else {
            index = (size_t)datamap_get(target->argid, longname, (int)strlen(longname), &nullflag);
            if (!nullflag && target->argvs[index]->argpresent)
                return argument_getf(target->argvs[index], argindex, defaultval);
        }
    }
    return defaultval;
}
char* args_getstr(args_t* target, const char* longname, int argindex, char* defaultval){
    size_t index;
    int nullflag=0;
    if (target) {
        if (!longname) {
            index = 0;
            return argument_gets(target->argvs[0], argindex, defaultval);
        }
        else {
            index = (size_t)datamap_get(target->argid, longname, (int)strlen(longname), &nullflag);
            if (!nullflag && target->argvs[index]->argpresent)
                return argument_gets(target->argvs[index], argindex, defaultval);
        }
    }
    return defaultval;
}
int args_getintfromstr(args_t* target, const char* longname, int argindex, int defaultval) {
    char* strres;
    int result;
    strres = args_getstr(target, longname, argindex, NULL);
    if (!strres) result = defaultval;
    else result = atoi(strres);
    return result;
}
double args_getdoublefromstr(args_t* target, const char* longname, int argindex, double defaultval) {
    char* strres;
    double result;
    strres = args_getstr(target, longname, argindex, NULL);
    if (!strres) result = defaultval;
    else result = atof(strres);
    return result;
}

int args_countargs(args_t* target, const char* longname) {
    size_t index;
    int nullflag = 0;
    if (target) {
        if (!longname) {
            index = 0;
            return target->argvs[index]->niargs + target->argvs[index]->nfargs + target->argvs[index]->nsargs;
        }
        else {
            index = (size_t)datamap_get(target->argid, longname, (int)strlen(longname), &nullflag);
            if (!nullflag && target->argvs[index]->argpresent)
                return target->argvs[index]->niargs + target->argvs[index]->nfargs + target->argvs[index]->nsargs;
        }
    }
    return 0;
}
int args_countfreeargs(args_t* target) {
    return args_countargs(target, NULL);
}
int args_ispresent(args_t* target, const char* longname) {
    size_t index;
    int nullflag;
    if (target) {
        index = (size_t)datamap_get(target->argid, longname, (int)strlen(longname), &nullflag);
        if (!nullflag) return target->argvs[index]->argpresent;
    }
    return 0;
}

void args_fprintgeneralhelp(FILE* target, args_t* arg, int colsperline) {
    size_t i,j;
    size_t maxargstrlen;
    size_t argstrlen;
    size_t remstr;
    size_t startpoint;
    size_t endpoint;
    size_t tmppoint;
    char* sdesc;
    if (target && arg) {
        maxargstrlen = 0;
        for (i = 1;i < arg->nargts; i++) {
            argstrlen = strlen(arg->argvs[i]->longname);
            if (argstrlen > maxargstrlen)
                maxargstrlen = argstrlen;
        }
        if (maxargstrlen > 20)maxargstrlen = 20;
        for (i = 1;i < arg->nargts; i++) {
            argstrlen = strlen(arg->argvs[i]->longname);
            if (arg->argvs[i]->shortname != 0)
                fprintf(target, "-%c or ", arg->argvs[i]->shortname);
            else
                fprintf(target, "        ");
            fprintf(target, "--%s  ", arg->argvs[i]->longname);
            while (argstrlen < maxargstrlen) {
                fputc(' ', target);
                argstrlen++;
            }
            startpoint = 0;
            if (arg->argvs[i]->shortdesc)
                sdesc = arg->argvs[i]->shortdesc;
            else
                sdesc = "< no short description available >";
            remstr = strlen(sdesc);
            endpoint = 0;
            while (endpoint < remstr) {
                /* find endpoint */
                tmppoint = startpoint;
                endpoint = startpoint;
                while (sdesc[tmppoint] != 0 && sdesc[tmppoint]!='\n' &&
                       tmppoint - startpoint < colsperline - maxargstrlen - 10)
                {
                    if (sdesc[tmppoint] == ' ')
                        endpoint = tmppoint;
                    tmppoint++;
                }
                if (sdesc[tmppoint] == ' ' || sdesc[tmppoint] == '\n' || sdesc[tmppoint] == 0)
                    endpoint = tmppoint;
                if (endpoint == startpoint && tmppoint > endpoint)endpoint = tmppoint - 1;
                /* print from startpoint to endpoint */
                for (j = startpoint;j <= endpoint;j++) {
                    fputc(sdesc[j], target);
                }
                /* add enough space for the next line to look good */
                if (endpoint < remstr) {
                    argstrlen = 0;
                    if(sdesc[tmppoint] != '\n')
                        fputc('\n', target);
                    while (argstrlen < maxargstrlen+10) {
                        fputc(' ', target);
                        argstrlen++;
                    }
                }
                startpoint = endpoint + 1;
            }
            fputc('\n', target);
        }
    }
}

void args_fprintspecifichelp(FILE* target, args_t* arg, const char* longname, int colsperline) {
    size_t argid;
    int nullflag;
    size_t i,j;
    size_t remstr;
    size_t endpoint;
    size_t startpoint;
    size_t tmppoint;
    char* fullname = "UNNAMED ARGUMENT";
    char* longdesc = "<Missing description>";
    if (target && arg && longname) {
        argid = (size_t)datamap_get(arg->argid, longname, (int)strlen(longname), &nullflag);
        if (!nullflag) {
            if (colsperline == -1 && args_ispresent(arg, "help"))
                colsperline = args_getint(arg,"help",0,80);
            if (arg->argvs[argid]->fullname)
                fullname = arg->argvs[argid]->fullname;
            if (arg->argvs[argid]->description)
                longdesc = arg->argvs[argid]->description;
            fprintf(target, "%s\n", fullname);
            /* --- */
            fprintf(target,"  Usage:\n    --%s", longname);
            i = 0;
            while (arg->argvs[argid]->argpattern[i] > 0) {
                if (arg->argvs[argid]->argpattern[i] == 'i')
                    fprintf(target, " <int>");
                if (arg->argvs[argid]->argpattern[i] == 'f')
                    fprintf(target, " <real>");
                if (arg->argvs[argid]->argpattern[i] == 's')
                    fprintf(target, " <text>");
                if (arg->argvs[argid]->argpattern[i] == 'v')
                    fprintf(target, " <variable input>");
                if (arg->argvs[argid]->argpattern[i] == '.')
                    fprintf(target, " <text> [<text>...]");
                i++;
            }
            if (arg->argvs[argid]->shortname != 0) {
                fprintf(target, "\n    -%c", arg->argvs[argid]->shortname);
                i = 0;
                while (arg->argvs[argid]->argpattern[i] > 0) {
                    if (arg->argvs[argid]->argpattern[i] == 'i')
                        fprintf(target, " <int>");
                    if (arg->argvs[argid]->argpattern[i] == 'f')
                        fprintf(target, " <real>");
                    if (arg->argvs[argid]->argpattern[i] == 's')
                        fprintf(target, " <text>");
                    if (arg->argvs[argid]->argpattern[i] == 'v')
                        fprintf(target, " <variable input>");
                    if (arg->argvs[argid]->argpattern[i] == '.')
                        fprintf(target, " <text> [<text>...]");
                    i++;
                }
            }
            /* --- */
            fprintf(target, "\n\n  Description:\n    ");
            remstr = strlen(longdesc);
            startpoint = 0;
            endpoint = 0;
            while (endpoint < remstr) {
                /* find endpoint */
                tmppoint = startpoint;
                endpoint = startpoint;
                while (longdesc[tmppoint] != 0 && longdesc[tmppoint] != '\n' &&
                    tmppoint - startpoint < (size_t)colsperline - 4)
                {
                    if (longdesc[tmppoint] == ' ')
                        endpoint = tmppoint;
                    tmppoint++;
                }
                if (longdesc[tmppoint] == ' ' || longdesc[tmppoint] == '\n' || longdesc[tmppoint] == 0)
                    endpoint = tmppoint;
                if (endpoint == startpoint && tmppoint>startpoint)endpoint = tmppoint - 1;
                /* print from startpoint to endpoint */
                for (j = startpoint;j <= endpoint;j++) {
                    fputc(longdesc[j], target);
                }
                /* add enough space for the next line to look good */
                if (endpoint < remstr) {
                    if (longdesc[tmppoint] != '\n')
                        fputc('\n', target);
                    fprintf(target, "    ");
                }
                startpoint = endpoint + 1;
            }
            fputc('\n', target);
        }
        else {
            fprintf(target, "[INVALID ARGUMENT]\n");
        }
    }
}

int args_get_quietness(args_t* target)
{
    return target->quiet_level;
}

void args_report_error(args_t* target, const char* msg, ...) {
    va_list va_p;
    static args_t* old_target = NULL;
    if (!target)target = old_target;
    else old_target = target;
    va_start(va_p, msg);
    if (!target || target->quiet_level < QUIETERR) {
        if (msg) {
            fprintf(stderr, "[ERROR] ");
            vfprintf(stderr, msg, va_p);
        }
    }
    va_end(va_p);
}
void args_report_warning(args_t* target, const char* msg, ...) {
    va_list va_p;
    static args_t* old_target = NULL;
    if (!target)target = old_target;
    else old_target = target;
    va_start(va_p, msg);
    if (!target || target->quiet_level < QUIETWARN) {
        if (msg) {
            fprintf(stderr, "[WARN] ");
            vfprintf(stderr, msg, va_p);
        }
    }
    va_end(va_p);
}
void args_report_progress(args_t* target, const char* msg, ...) {
    va_list va_p;
    static args_t* old_target = NULL;
    if (!target)target = old_target;
    else old_target = target;
    va_start(va_p, msg);
    if (!target || target->quiet_level < QUIETPROGRESS) {
        if (msg) {
            vfprintf(stderr, msg, va_p);
        }
    }
    va_end(va_p);
}
void args_report_info(args_t* target, const char* msg, ...) {
    va_list va_p;
    static args_t* old_target = NULL;
    if (!target)target = old_target;
    else old_target = target;
    va_start(va_p, msg);
    if (!target || target->quiet_level < QUIETINFO) {
        if (msg) {
            fprintf(stderr, "[INFO] ");
            vfprintf(stderr, msg, va_p);
        }
    }
    va_end(va_p);
}

int args_is_helpmode(args_t* target) {
    return target->helpmodeon;
}