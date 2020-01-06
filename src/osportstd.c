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
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "textparsing.h"
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
char* os_extractdirname(char* path) {
    char* result;
    size_t i;
    i = strlen(path);
    while (i > 0 && path[i - 1] != '/' && path[i - 1] != '\\') {
        i--;
    }
    if (i == 0) result = "";
    else {
        result =(char*)  malloc(i + 1);
        memcpy(result, path, i);
        result[i] = 0;
    }
    return result;
}
#ifdef _WIN32
/* helper functions for windows-specific prcoess output parsing */
size_t _win32_ReadFromPipe(HANDLE g_hChildStd_OUT_Rd, char* result, size_t buffersize, int* nullflag, int wait) {
    DWORD readbytes, bytesavailable, bltm;
    char tmpc;
    int tmp;
    readbytes = 0;
    OVERLAPPED sto = { 0 };

    if (!wait)
        PeekNamedPipe(g_hChildStd_OUT_Rd, &tmpc, 1, &readbytes, &bytesavailable, &bltm);
    else
        bytesavailable = 1;
    if (bytesavailable == 0 || !ReadFile(g_hChildStd_OUT_Rd, result, (int32_t)buffersize - 1, &readbytes, &sto)) {
        *nullflag = 1;
    }
    else {
        *nullflag = 0;
    }
    tmp = GetLastError();
    result[readbytes] = 0;
    return readbytes;
}
HANDLE _win32_CreateChildProcess(char* prog, int argc, char** argv,HANDLE g_hChildStd_IN_Rd, HANDLE g_hChildStd_IN_Wr, HANDLE g_hChildStd_OUT_Rd, HANDLE g_hChildStd_OUT_Wr) {
    PROCESS_INFORMATION procinfo;
    STARTUPINFO startinfo;
    SECURITY_ATTRIBUTES secattr;
    HANDLE hresult;
    char eof = EOF;

    hresult = INVALID_HANDLE_VALUE;

    char* cmdline;
    cmdline = text_join(argv, "\" \"", "\"","\"", argc);
    /* note : for the sake of clarity I avoid using microsoft macros when a serviceable ANSI-C alternative exists */
    memset(&procinfo, 0, sizeof(PROCESS_INFORMATION));
    memset(&startinfo, 0, sizeof(STARTUPINFO));
    startinfo.cb = sizeof(STARTUPINFO);
    startinfo.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    startinfo.hStdOutput = g_hChildStd_OUT_Wr;
    startinfo.hStdInput = g_hChildStd_IN_Rd;
    startinfo.dwFlags |= STARTF_USESTDHANDLES;
    memset(&secattr, 0, sizeof(SECURITY_ATTRIBUTES));
    secattr.bInheritHandle = 1;
    if (!CreateProcess(prog, cmdline, &secattr, &secattr, 1, 0, NULL, NULL, &startinfo, &procinfo)) {
        fprintf(stderr, "Process did not succeed: %s\n", cmdline);
    }
    else {
        CloseHandle(g_hChildStd_OUT_Wr);
        CloseHandle(g_hChildStd_IN_Rd);
        CloseHandle(procinfo.hThread);
        WaitForInputIdle(procinfo.hProcess, INFINITE);
        hresult = procinfo.hProcess;
    }
    return hresult;
}
#endif

char* os_stdoutfromexec(char* program, int argc, char** argv, size_t* outlen, size_t bufsize) {
    char* buffer;
    char* result;
    size_t totlen;
    size_t curlen;
    size_t oldlen;
#ifdef _WIN32
    SECURITY_ATTRIBUTES sattr;
    int nf;

    HANDLE g_hChildStd_IN_Rd = NULL;
    HANDLE g_hChildStd_IN_Wr = NULL;
    HANDLE g_hChildStd_OUT_Rd = NULL;
    HANDLE g_hChildStd_OUT_Wr = NULL;
    HANDLE hprochandle;
    int32_t exitstatus;

    sattr.nLength = sizeof(SECURITY_ATTRIBUTES);
    sattr.bInheritHandle = TRUE;
    sattr.lpSecurityDescriptor = NULL;
    /*
    Initialize the pipes and make sure the child does not inherit write/read access
    to the wrong end of each pipe. Then launch the external process.
    */
    if (!CreatePipe(&g_hChildStd_OUT_Rd, &g_hChildStd_OUT_Wr, &sattr, 0))
        fprintf(stderr, "Could not open Pipe to child process STDOUT\n");
    if (!SetHandleInformation(g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0))
        fprintf(stderr, "Could not prevent inheritance of child STDOUT read access\n");
    if (!CreatePipe(&g_hChildStd_IN_Rd, &g_hChildStd_IN_Wr, &sattr, 0))
        fprintf(stderr, "Could not open Pipe to child process STDIN\n");
    if (!SetHandleInformation(g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0))
        fprintf(stderr, "Could not prevent inheritance of child STDIN write access\n");
    hprochandle = _win32_CreateChildProcess(program, argc, argv, g_hChildStd_IN_Rd, g_hChildStd_IN_Wr, g_hChildStd_OUT_Rd, g_hChildStd_OUT_Wr);

    nf = 0;
    totlen = oldlen = 0;
    buffer = (char*)calloc(bufsize, 1);
    curlen = _win32_ReadFromPipe(g_hChildStd_OUT_Rd, buffer, bufsize, &nf, 0);
    result = (char*)malloc(curlen+1);
    GetExitCodeProcess(hprochandle, &exitstatus);
    /*
    nf can be raised if the process does not generate its stdout output quickly enough.
    Also, some programs will not termina before their ouput has been parsed.
    This explains the conditions below
    */
    while (!nf || exitstatus == STILL_ACTIVE) {
        curlen = _win32_ReadFromPipe(g_hChildStd_OUT_Rd, buffer, bufsize, &nf, 0);
        if (curlen > 0) {
            totlen += curlen;
            result = (char*)realloc(result, totlen + 1);
            memcpy(result + oldlen, buffer, curlen);
            oldlen = totlen;
        }
        GetExitCodeProcess(hprochandle, &exitstatus);
        if (nf && exitstatus != STILL_ACTIVE)Sleep(1);
    }
    curlen = _win32_ReadFromPipe(g_hChildStd_OUT_Rd, buffer, bufsize, &nf, 0);
    if (curlen > 0) {
        totlen += curlen;
        result = (char*)realloc(result, totlen + 1);
        memcpy(result + oldlen, buffer, curlen);
    }
    free(buffer);
    result[totlen] = 0;
#else
    int exitstatus;
    int pipefd[6];
    pid_t cpid;
    char** nargv;
    result = NULL;
    if (argc < 1) return NULL;

    if (pipe(pipefd) == -1) {
        fprintf(stderr, "could not create stdout pipe to %s:\n", program);
        return NULL;
    }
    if (pipe(pipefd+2) == -1) {
	close(pipefd[0]); close(pipefd[1]);
        fprintf(stderr, "could not create stderr pipe to %s\n", program);
        return NULL;
    }
    if (pipe(pipefd+4) == -1) {
	close(pipefd[0]); close(pipefd[1]);
	close(pipefd[2]); close(pipefd[3]);
        fprintf(stderr, "could not create stdin pipe to %s\n", program);
        return NULL;
    }
    /* write pipe is not used */
    cpid = fork();
    if (cpid == -1) {
        fprintf(stderr, "could not fork process\n");
    }
    else if (cpid != 0) {
        /* parent process */
        close(pipefd[1]);
        close(pipefd[3]);
        close(pipefd[4]);

        totlen = oldlen = 0;
        buffer = (char*)calloc(bufsize, 1);
        curlen = 1000;
        curlen = read(pipefd[0], buffer, bufsize - 1);
        while ( curlen > 0) {
            totlen += curlen;
            result = (char*)realloc(result, totlen + 1);
            memcpy(result + oldlen, buffer, curlen);
            result[totlen] = 0;
            oldlen = totlen;
            curlen = read(pipefd[0], buffer, bufsize - 1);
            if(waitpid(cpid, &exitstatus, WNOHANG))break;
        }
        while ( curlen > 0) {
            totlen += curlen;
            result = (char*)realloc(result, totlen + 1);
            memcpy(result + oldlen, buffer, curlen);
            result[totlen] = 0;
            oldlen = totlen;
            curlen = read(pipefd[0], buffer, bufsize - 1);
        }
        if(!result || strlen(result)<1) {
            fprintf(stderr,"program produced no output\n");
        }
        free(buffer);
        close(pipefd[0]);
        close(pipefd[2]);
	close(pipefd[5]);
        wait(0);
    }
    else {
        /* child process */
        close(pipefd[0]);
	close(pipefd[2]);
        close(pipefd[5]);
        while ((dup2(pipefd[1], STDOUT_FILENO) == -1) && (errno == EINTR)) {}
        while ((dup2(pipefd[3], STDERR_FILENO) == -1) && (errno == EINTR)) {}
        while ((dup2(pipefd[4], STDIN_FILENO) == -1) && (errno == EINTR)) {}

        nargv=malloc(sizeof(char*)*(argc+1));
        memcpy(nargv,argv,sizeof(char*)*(argc));
        nargv[argc]=NULL;
        exitstatus = execv(nargv[0], nargv);
	if(exitstatus){
            fprintf(stdout,"program exited with status [%d], Arguments (%d) were:\n",errno, argc);
            curlen = 0;
/*
            while(nargv[curlen] != NULL) {
                fprintf(stdout,"%s ",nargv[curlen++]);
            }
            fprintf(stdout,"\n");*/
        }
        close(pipefd[1]);
        close(pipefd[3]);
	close(pipefd[4]);
        free(nargv);
        _exit((char)errno);
    }
#endif
    *outlen = totlen;
    return result;
}

void* memcpyalloc(void* data, size_t datasize) {
    void* result;
    result = malloc(datasize);
    memcpy(result, data, datasize);
    return result;
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
            copysize = f->maxbpos - (size_t)(f->bufpos);
        memcpy(cbuf + curstart, f->buffer + f->bufpos, copysize);
        curstart += copysize;
        f->bufpos += copysize;
    }
    while (curstart < resultsize && f->bufpos >= (long long)f->maxbpos && f->maxbpos == (long long)PFBUFFER - 1) {
#ifdef _WIN32
        f->maxbpos = _fread_nolock(f->buffer, 1, PFBUFFER - 1, f->filepointer);
        tmp = ferror(f->filepointer);
#else
        f->maxbpos = fread(f->buffer, 1, PFBUFFER - 1, f->filepointer);
#endif
        f->bufpos = 0;
        f->buffer[f->maxbpos] = 0;
        /* we assume that the previous string ends with a \0 and therefore move the caret back */
        if (resultsize - curstart > f->maxbpos)
            copysize = f->maxbpos;
        else
            copysize = resultsize - curstart;
        memcpy(cbuf + curstart, f->buffer + f->bufpos, copysize);
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
    return _ftelli64(f->filepointer) + f->bufpos - PFBUFFER;
#else
    return (long long)ftell(f->filepointer) + f->bufpos - PFBUFFER;
#endif
#else
    return (long long)ftell(f->filepointer) + f->bufpos - PFBUFFER;
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
    return fseek(f->filepointer, (long)offset, origin);
#endif
#else
    return fseek(f->filepointer, offset, origin);
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
int PFseemslikeASCIItext(PF_t* f, size_t bytes_to_test) {
    char* buffer;
    long long start;
    size_t i;
    int result;
    if (!f || !(f->filepointer)) return 2;
    buffer = (char*)malloc(bytes_to_test);
    #ifdef _WIN32
    #ifdef _WIN64
    start = _ftelli64(f->filepointer);
    rewind(f->filepointer);
    bytes_to_test = _fread_nolock(buffer, 1, bytes_to_test, f->filepointer);
    _fseeki64(f->filepointer, start, SEEK_SET);
    #else
    start = (long long)ftell(f->filepointer);
    rewind(f->filepointer);
    bytes_to_test = _fread_nolock(buffer, 1, bytes_to_test, f->filepointer);
    fseek(f->filepointer, (long)start, SEEK_SET);
    #endif
    #else
    start = (long long)ftell(f->filepointer);
    rewind(f->filepointer);
    bytes_to_test = fread(buffer, 1, bytes_to_test, f->filepointer);
    fseek(f->filepointer, start, SEEK_SET);
    #endif
    result = 1;
    for (i = 0;i < bytes_to_test;i++) {
        if (buffer[i] == '\t' || buffer[i] == '\n' || buffer[i] == '\r') continue;
        if (buffer[i] < ' ') {
            result = 0;
            break;
        }
        if (result > '~') {
            result = 0;
            break;
        }
    }
    free(buffer);
    return result;
}
char* PFgetfullpathptr(PF_t* f) {
    return f->filename;
}
char* PFgetbasenameptr(PF_t* f) {
    return f->filename + f->basenamestart;
}
FILE* PFgetFILE(PF_t* f) {
    if (!f)return NULL;
    return f->filepointer;
}
size_t PFnumlines(PF_t* f, int countempty) {
    ssize_t start;
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
    start=(ssize_t)ftell(f->filepointer);
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
    bufpos = (size_t) f->bufpos;
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
unimem PFreadunimem16(PF_t* f) {
    uint16_t tmpsz;
    unimem output;
    tmpsz = PFgetint16(f);
    output = (unimem)malloc((size_t)tmpsz);
    PFread(((char*)output) + 1, 1, tmpsz, f);
    return output;
}
unimem PFreadunimem32(PF_t* f) {
    uint32_t tmpsz;
    unimem output;
    tmpsz = PFgetint32(f);
    output = (unimem)malloc((size_t)tmpsz);
    PFread(((char*)output) + 1, 1, tmpsz, f);
    return output;
}
unimem PFreadunimem64(PF_t* f) {
    uint64_t tmpsz;
    unimem output;
    tmpsz = PFgetint64(f);
    output = (unimem)malloc((size_t)tmpsz);
    PFread(((char*)output) + 1, 1, (size_t)tmpsz, f);
    return output;
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