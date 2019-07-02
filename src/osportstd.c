#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "osportstd.h"

int os_fileexists(const char* path) {
#ifdef _WIN32
    if (_access(path, 0) != -1)
        return 1;
#else
    if (access(path, F_OK) != -1)
        return 1;
#endif
    return 0;
}
char* os_rmdirname(char* path) {
    size_t i;
    i = strlen(path);
    while (i > 0 && path[i - 1] != '/' && path[i - 1] != '\\') {
        i--;
    }
    return path + i;
}

struct portablefile_t {
    FILE* filepointer;
    char* filename;
    int dofclose;
    size_t fnlen;
    size_t basenamestart;
    size_t bnlen;
    char buffer[PFBUFFER];
    size_t maxbpos;
    long long bufpos;
};

/* adapted std calls*/
errcode_t PFopen(PF_t** f, const char* path, const char* mode) {
    errcode_t outval;
    PF_t* result;
    int notfound;
    size_t cpos;
    outval = 0;
    result = calloc(1,sizeof(PF_t));
    result->dofclose = 1;
    if (!result) outval = 1;
    else {
#if PFSTDIOENABLE == 1
        if (strcmp(path, "stdout") == 0) {
            result->filepointer = stdout;
            result->dofclose = 0;
        }
        else if (strcmp(path, "stdin") == 0){
            result->filepointer = stdin;
            result->dofclose = 0;
        }
        else if (strcmp(path, "stderr") == 0){
            result->filepointer = stderr;
            result->dofclose = 0;
        }
        else 
#endif
#ifdef _WIN32
        outval = (errcode_t) fopen_s(&(result->filepointer), path, mode);
#else
        result->filepointer = fopen(path, mode);
        if (!result->filepointer) outval = 2;
#endif
        result->fnlen = strlen(path);
        result->filename = malloc(result->fnlen + 1);
        memcpy(result->filename, path, result->fnlen);
        result->filename[result->fnlen] = 0;
        result->bnlen = 0;
        notfound = 1;
        cpos = result->fnlen;
        while (cpos > 0 && notfound) {
            cpos--;
#ifdef _WIN32
            if (result->filename[cpos] == '/' || result->filename[cpos] == '\\')
                notfound = 0;
#else
            if (result->filename[cpos] == '/')
                notfound = 0;
#endif
        }
        if (notfound)
            result->basenamestart = 0;
        else
            result->basenamestart = cpos+1;
        result->bnlen = result->fnlen - result->basenamestart;
        result->bufpos = PFBUFFER;
        result->maxbpos = PFBUFFER-1;
    }
    *f = result;
    return outval;
}
errcode_t PFclose(PF_t* f) {
    if (f) {
        if (f->filepointer && f->dofclose)fclose(f->filepointer);
        if (f->filename)free(f->filename);
        memset(f, 0, sizeof(PF_t));
        free(f);
        return 0;
    }
    return 1;
}
size_t PFread(void* buffer, size_t element_size, size_t element_count, PF_t* f) {
    size_t copysize;
    size_t resultsize;
    size_t curstart;
#ifdef _WIN32
    int tmp;
#endif
    char* cbuf;
    if (!f || !(f->filepointer)) {
        return 0;
    }
    cbuf = (char*)buffer;
    resultsize = element_size*element_count;
    curstart = 0;
    copysize = resultsize;
    if (f->bufpos < (long long)f->maxbpos) {
        if (f->bufpos + copysize > f->maxbpos)
            copysize = f->maxbpos - f->bufpos + 1;
        memcpy(cbuf + curstart, f->buffer + f->bufpos, copysize);
        if (copysize == f->maxbpos + 1)
            copysize--;
        curstart += copysize;
        f->bufpos += copysize;
    }
    while (f->bufpos >= (long long)f->maxbpos && f->maxbpos == (long long)PFBUFFER - 1) {
#ifdef _WIN32
        f->maxbpos = _fread_nolock(f->buffer, 1, PFBUFFER - 1, f->filepointer);
        tmp = ferror(f->filepointer);
#else
        f->maxbpos = fread(f->buffer, 1, PFBUFFER - 1, f->filepointer);
#endif
        f->bufpos = 0;
        f->buffer[f->maxbpos] = 0;
        /* we assume that the previous string ends with a \0 and therefore move the caret back */
        if (curstart > 0) curstart--;
        if (resultsize - curstart > f->maxbpos)
            copysize = f->maxbpos + 1;
        else
            copysize = resultsize - curstart;
        memcpy(cbuf + curstart, f->buffer + f->bufpos, copysize);
        if (copysize == f->maxbpos + 1)
            copysize--;
        curstart += copysize;
        f->bufpos += copysize;
    }
    return curstart;
}
size_t PFwrite(const void* buffer, size_t element_size, size_t element_count, PF_t* f) {
    if (f && f->filepointer) {
        return fwrite(buffer, element_size, element_count, f->filepointer);
    }
    return 0;
}
long long PFtell(PF_t* f) {
    if (!f || !(f->filepointer)) return -1;
#ifdef _WIN32
#ifdef _WIN64
    return _ftelli64(f->filepointer) + f->bufpos;
#else
    return (long long)ftell(f->filepointer) + f->bufpos;
#endif
#else
    return (long long)ftell(f->filepointer) + f->bufpos;
#endif
}
errcode_t PFseek(PF_t* f, long long offset, int origin) {
    if (!f)return 1;
    f->maxbpos = PFBUFFER - 1;
    f->bufpos = PFBUFFER;
#ifdef _WIN32
#ifdef _WIN64
    return _fseeki64(f->filepointer,offset, origin);
#else
    return (long long)ftell(f->filepointer) + f->bufpos;
#endif
#else
    return (long long)ftell(f->filepointer) + f->bufpos;
#endif
}
errcode_t PFrewind(PF_t* f) {
    if (!f)return 1;
    f->maxbpos = PFBUFFER - 1;
    f->bufpos = PFBUFFER;
    rewind(f->filepointer);
    return 0;
}
int PFgetc(PF_t* f) {
    if (!f)return -1;
    if (f->bufpos < (long long)f->maxbpos) {
        f->bufpos++;
        return f->buffer[f->bufpos - 1];
    }
    if (f->maxbpos < PFBUFFER - 1) return -1;
#ifdef _WIN32
    f->maxbpos = fread_s(f->buffer, PFBUFFER, 1, PFBUFFER - 1, f->filepointer);
#else
    f->maxbpos = fread(f->buffer, 1, PFBUFFER - 1, f->filepointer);
#endif
    if (f->maxbpos == 0)return -1;
    f->bufpos = 1;
    return f->buffer[0];
}
errcode_t PFputc(int c, PF_t* f) {
    if (!f)return -1;
    return fputc(c,f->filepointer);
}
errcode_t PFprintf(PF_t* f, const char* fmt, ...) {
    errcode_t res;
    va_list argp;
    va_start(argp, fmt);
    if (!f || !(f->filepointer))return -1;
    res = (errcode_t) vfprintf(f->filepointer, fmt, argp);
    va_end(argp);
    return res;
}
errcode_t PFvprintf(PF_t* f, const char* fmt, va_list argp) {
    return vfprintf(f->filepointer, fmt, argp);
}
/* extra functions */
FILE* PFgetFILE(PF_t* f) {
    if (!f)return NULL;
    return f->filepointer;
}
size_t PFnumlines(PF_t* f, int countempty) {
    long long start;
    size_t nread;
    size_t i;
    size_t nlines;
    int stop;
    int nonempty;
    char buffer[0x10000] = { 0 };
    if (!f || !f->filepointer) return 0;
#ifdef _WIN64
    start=_ftelli64(f->filepointer);
#else
    start=(long long)ftell(f->filepointer);
#endif
    rewind(f->filepointer);
    stop = 0;
    nonempty = 0;
    nlines = 0;
    while (!stop) {
        nread = fread(buffer, 1, 0x10000-1, f->filepointer);
        buffer[nread] = 0;
        for (i = 0;i < nread;i++) {
            if (buffer[i] == '\n') {
                if (countempty || nonempty)nlines++;
                nonempty = 0;
            }
            else if (nonempty == 0 && buffer[i] != '\r') {
                nonempty = 1;
            }
        }
        if (nread < 0x10000-1)stop = 1;
    }
    if (nonempty)nlines++;
#ifdef _WIN64
    _fseeki64(f->filepointer, start, SEEK_SET);
#else
    fseek(f->filepointer, start, SEEK_SET);
#endif
    return nlines;
}

int16_t PFgetint16(PF_t* f) {
    int16_t res=0;
    if(f && f->filepointer)
        PFread(&res, sizeof(int16_t), 1, f);
    return res;
}
int32_t PFgetint32(PF_t* f) {
    int32_t res = 0;
    if (f && f->filepointer)
        PFread(&res, sizeof(int32_t), 1, f);
    return res;
}
int64_t PFgetint64(PF_t* f) {
    int64_t res = 0;
    if (f && f->filepointer)
        PFread(&res, sizeof(int64_t), 1, f);
    return res;
}
float PFgetfloat(PF_t* f) {
    float res = 0.0;
    if (f && f->filepointer)
        PFread(&res, sizeof(float), 1, f);
    return res;
}
double PFgetdouble(PF_t* f) {
    double res = 0.0;
    if (f && f->filepointer)
        PFread(&res, sizeof(double), 1, f);
    return res;
}
char* _PFstrreadline(const char* buffer, size_t* start, size_t* extension) {
    size_t end;
    size_t len;
    char* result;
    end = *start;
    while (buffer[end] != 0 && buffer[end] != '\n') {
        end++;
    }
    len = end - *start;
    if (len > 0 && buffer[end - 1] == '\r')len--;
    result = malloc(len + 1);
    memcpy(result, buffer + *start, len);
    result[len] = 0;
    *start = end + 1;
    if (extension) *extension = len;
    return result;
}
char* PFreadline(PF_t* f) {
    size_t extension;
    size_t totlinelen;
    size_t readbytes;
    size_t bufpos;
    char* line;
    char* lineend;
    if (!f || !f->filepointer) return NULL;
    if (f->bufpos < PFBUFFER && f->buffer[f->bufpos] == 0)return NULL;
    line = calloc(1, sizeof(char));
    bufpos = f->bufpos;
    totlinelen = 0;
    if (bufpos <= PFBUFFER-1) {
        lineend = _PFstrreadline(f->buffer, &bufpos, &extension);
        totlinelen = extension;
        line = realloc(line, totlinelen + 1);
        memcpy(line, lineend, extension);
        free(lineend);
    }
    while (bufpos>PFBUFFER-1) {
        bufpos = 0;
#ifdef _WIN32
        readbytes = fread_s(f->buffer, PFBUFFER, 1, PFBUFFER-1, f->filepointer);
#else
        readbytes = fread(f->buffer, 1, PFBUFFER-1, f->filepointer);
#endif
        f->buffer[readbytes] = 0;
        lineend = _PFstrreadline(f->buffer, &bufpos, &extension);
        totlinelen += extension;
        line = realloc(line, totlinelen + 1);
        memcpy(line + totlinelen - extension, lineend, extension);
        free(lineend);
    }
    if (bufpos > 0 && f->buffer[bufpos - 1] == 0)
        f->bufpos = bufpos - 1;
    else
        f->bufpos = bufpos;
    line[totlinelen] = 0;
    return line;
}

errcode_t PFputint16(PF_t* f, int16_t input){
    if (PFwrite(&input, sizeof(int16_t), 1, f) != sizeof(int16_t))return 1;
    return 0;
}
errcode_t PFputint32(PF_t* f, int32_t input){
    if (PFwrite(&input, sizeof(int32_t), 1, f) != sizeof(int32_t))return 1;
    return 0;
}
errcode_t PFputint64(PF_t* f, int64_t input){
    if (PFwrite(&input, sizeof(int64_t), 1, f) != sizeof(int64_t))return 1;
    return 0;
}
errcode_t PFputfloat(PF_t* f, float input){
    if (PFwrite(&input, sizeof(float), 1, f) != sizeof(float))return 1;
    return 0;
}
errcode_t PFputdouble(PF_t* f, double input){
    if (PFwrite(&input, sizeof(double), 1, f) != sizeof(double))return 1;
    return 0;
}