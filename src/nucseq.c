#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nucseq.h"
#include "vecops.h"
#include "hash_table.h"
#include "datalist.h"

inline int _nucmatch(nucleotide from, nucleotide to) {
    // if from = nucN, any nucleotide returns a match
    // if to = nucN, we don't really know
    // return values
    // 0: FALSE, 1: TRUE, 2:MAYBE
    if (from == nucN) return 1;
    if (from == to) return 1;
    if (to == nucN) return 2;
    return 0;
}

inline nucleotide _getnuc(nucseq* seq, long long pos) {
    if (pos >= 0 && pos < (long long)seq->len)
        return seq->seq[pos];
    else
        return nucN;
}

nucleotide nuc_complement(nucleotide nuc)
{
    if (nuc >= 0)
        return 3 - nuc;
    else
        return -1;
}

const int c2ntable[] = {
    //        0  1  2  3  4  5  6  7  8  9  A  B  C  D  E  F
    /*0*/    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    /*1*/    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    /*2*/    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    /*3*/    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    /*4*/    -1, nucA,-1, nucC,-1,-1,-1, nucG,-1,-1,-1,-1,-1,-1,-1,-1,
    /*5*/    -1,-1,-1,-1, nucT, nucT,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    /*6*/    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    /*7*/    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
const char n2ctable[] = { 'A', 'C', 'G', 'T' };

nucleotide char2nuc(char c) {
    return c2ntable[c];
}

char nuc2char(nucleotide c) {
    return (c < 0) ? 'N' : n2ctable[c];
}

char trinuc2aa[] = {
    'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I',
    'Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L',
    'E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
    'X','Y','X','Y','S','S','S','S','X','C','W','C','L','F','L','F'
};

char translate_nuc(nucleotide first, nucleotide second, nucleotide wobble) {
    if (first > 3 || second > 3 || wobble > 3) return 'X';
    else return trinuc2aa[(first << 4) + (second << 2) + wobble];
}

void clear_nucseq(nucseq* target)
{
    if (!(target->flags&NOALLOC))
    {
        if (target->seq) free(target->seq);
        if (target->q) free(target->q);
    }
    target->seq = NULL;
    target->q = NULL;
    target->len = 0;
    target->flags = 0;
    if (target->name)free(target->name);
    target->name = NULL;
}
void force_clear_nucseq(nucseq* target)
{
    target->seq = NULL;
    target->q = NULL;
    target->len = 0;
    target->flags = 0;
    target->name = NULL;
}

size_t nucseq_from_string(nucseq* target, char* str)
{
    size_t i_;
    size_t reslen = 0;
    if (!target)return 0;
    if (!(target->flags&READONLY)) {
        clear_nucseq(target);
        // Read the string and allocate length
        while (str[reslen] != 0)
            reslen++;
        target->seq = (char*)malloc(sizeof(char)*reslen);
        if (!target->seq) {
            fprintf(stderr, "Not enough memory for sequence");
        }
        memcpy(target->seq, str, reslen);
        for (i_ = 0;i_ < reslen;i_++)
            target->seq[i_] = c2ntable[target->seq[i_]];
        target->len = reslen;
    }
    return target->len;
}
char* nucseq_tritext(nucseq* target) {
    char* result;
    size_t i;

    result = malloc(target->len - 1);
    result[0] = target->seq[0] * 16 + target->seq[1] * 4 + target->seq[2];
    for (i = 3; i < target->len; i++) {
        result[i - 2] = ((result[i - 3] << 2) & 0b00111111) + target->seq[i];
    }
    for (i = 0;i < target->len - 2;i++) {
        result[i] += ' ';
    }
    result[i] = 0;
    return result;
    return result;
}

char* strreadline(const char* buffer, size_t* start, size_t* extension) {
    size_t end;
    size_t len;
    char* result;
    end = *start;
    while (buffer[end] != 0 && buffer[end] != '\n') {
        end++;
    }
    len = end - *start;
    result = malloc(len + 1);
    memcpy(result, buffer + *start, len);
    result[len] = 0;
    *start = end + 1;
    if (extension) *extension = len;
    return result;
}

nucseq** nucseq_array_from_fasta(PF_t* f, size_t* OUT_count, int saveseqnames, size_t minlen) {
    nucseq** output = NULL;
    size_t nseqs;
    size_t nalloc;
    char* line;
    size_t stralloc = 0;
    size_t strfill = 0;
    size_t strcur = 0;
    size_t bufpos;
    char* curseqname = NULL;
    char* curseq = NULL;
    char buffer[0x10000] = { 0 };
    nalloc = 1;
    output = (nucseq**)malloc(sizeof(nucseq*));

    bufpos = 0x10000;
    line = PFreadline(f);
    nseqs = 0;
    while (line) {
        strcur = strlen(line);
        if (line[0] == '>') {
            if (curseq) {
                if (strfill >= minlen) {
                    nseqs++;
                    if (nalloc <= nseqs) {
                        nalloc *= 2;
                        output = (nucseq**)realloc(output, sizeof(nucseq*)*nalloc);
                    }
                    output[nseqs] = NULL;
                    output[nseqs - 1] = (nucseq*)malloc(sizeof(nucseq));
                    force_clear_nucseq(output[nseqs - 1]);
                    nucseq_from_string(output[nseqs - 1], curseq);
                    if (saveseqnames) {
                        output[nseqs - 1]->name = curseqname;
                    }
                }
                free(curseq);
                curseq = NULL;
                strfill = 0;
                stralloc = 0;
            }
            if (saveseqnames) {
                curseqname = malloc(strcur + 1);
                memcpy(curseqname, line + 1, strcur - 1);
                curseqname[strcur] = 0;
            }
            free(line);
        }
        else {
            if (strfill + strcur > stralloc) {
                stralloc = (strfill + strcur) * 2;
                curseq = (char*)realloc(curseq, stralloc + 1);
                curseq[stralloc] = 0;
            }
            memcpy(curseq + strfill, line, strcur);
            strfill += strcur;
            curseq[strfill] = 0;
            free(line);
        }
        line = PFreadline(f);
    }
    if (curseq) {
        if (strfill >= minlen) {
            nseqs++;
            if (nalloc <= nseqs) {
                nalloc *= 2;
                output = (nucseq**)realloc(output, sizeof(nucseq*)*nalloc);
            }
            output[nseqs] = NULL;
            output[nseqs - 1] = (nucseq*)malloc(sizeof(nucseq));
            force_clear_nucseq(output[nseqs - 1]);
            nucseq_from_string(output[nseqs - 1], curseq);
            if (saveseqnames) {
                output[nseqs - 1]->name = curseqname;
            }
        }
        free(curseq);
        curseq = NULL;
        strfill = 0;
        stralloc = 0;
    }
    output = (nucseq**)realloc(output, sizeof(nucseq*)*(nseqs + 1));
    *OUT_count = nseqs;
    return output;
}

void nucseq_array_clear(nucseq** target, size_t count) {
    size_t i;
    for (i = 0;i < count;i++) {
        if (target[i]) {
            clear_nucseq(target[i]);
            free(target[i]);
            target[i] = NULL;
        }
    }
    free(target);
}
char* nucseq2fasta(nucseq* target, const char* name)
{
    char* result = NULL;
    size_t namelen;
    size_t i_, curpos;
    namelen = strlen(name);
    result = (char*)malloc(1 + namelen + 1 + target->len + 2);
    if (result) {
        result[0] = '>';
        memcpy(result + 1, name, namelen);
        result[namelen + 1] = '\n';
        curpos = namelen + 2;
        for (i_ = 0;i_ < target->len;i_++)
            result[curpos++] = nuc2char(target->seq[i_]);
        result[curpos] = '\n';
        result[curpos + 1] = '\0';
    }
    else {
        fprintf(stderr, "Allocation of sequence failed");
    }
    return result;
}
void nucseq2fastafile(nucseq* target, PF_t* f, const char* name) {
    size_t i_;
    static int32_t altid = 0;
    if (name && strcmp(name, "0") == 0) {
        altid = 0;
    }
    else {
        if (name) {
            PFputc('>', f);
            PFwrite(name, 1, strlen(name), f);
        }
        else if (target->name) {
            PFputc('>', f);
            PFwrite(target->name, 1, strlen(target->name), f);
        }
        else {
            PFprintf(f, ">%d", ++altid);
        }
        PFputc('\n', f);
        for (i_ = 0;i_ < target->len;i_++)
            PFputc(nuc2char(target->seq[i_]), f);
        PFputc('\n', f);
    }
}
void nucseq_array2fastafile(nucseq** nucseq_array, size_t nseq, PF_t* f, int reset_seqid) {
    size_t i;
    if (reset_seqid) nucseq2fastafile(NULL, NULL, "0");
    for (i = 0;i < nseq; i++) {
        nucseq2fastafile(nucseq_array[i], f, NULL);
    }
}

void nucseq2fastqfile(nucseq* target, PF_t* f, const char* name) {
    size_t i_;
    PFputc('@', f);
    PFwrite(name, 1, strlen(name), f);
    PFputc('\n', f);
    for (i_ = 0;i_ < target->len;i_++)
        PFputc(nuc2char(target->seq[i_]), f);
    PFprintf(f,"\n+\n");
    if (target->q) {
        for (i_ = 0;i_ < target->len;i_++)
            PFputc(nuc2char('!' + target->q[i_]), f);
    }
    else {
        for (i_ = 0;i_ < target->len;i_++)
            PFputc('!' + 20, f);
    }
    PFputc('\n', f);
}
void nucseqrevcomp(nucseq * src, nucseq * target)
{
    size_t i_;
    if (target && !(target->flags & READONLY) && src)
    {
        clear_nucseq(target);
        target->seq = (char*)malloc(src->len);
        target->len = src->len;
        for (i_ = 0; i_ < src->len; i_++) {
            target->seq[target->len - i_ - 1] = nuc_complement(src->seq[i_]);
        }
        if (src->q) {
            target->q = (char*)malloc(src->len);
            for (i_ = 0; i_ < src->len; i_++) {
                target->q[target->len - i_ - 1] = src->q[i_];
            }
        }
        else target->q = NULL;
    }
}
void subseq(nucseq* target, nucseq* src, size_t start, size_t len)
{
    if (target && !(target->flags & READONLY) && src)
    {
        clear_nucseq(target);
        target->seq = (char*)malloc(len);
        if (!target->seq) {
            fprintf(stderr, "Not enough memory for sub-sequence");
        }
        memcpy(target->seq, src->seq + start, len);
        if (src->q)
        {
            target->q = (char*)malloc(len);
            if (!target->seq) {
                fprintf(stderr, "Not enough memory for sub-sequence q");
            }
            memcpy(target->q, src->q + start, len);
        }
        target->len = len;
    }
}
void subseq_nocpy(nucseq* target, nucseq* src, size_t start, size_t len)
{
    if (target && !(target->flags & READONLY) && src)
    {
        //clear_nucseq(target);
        target->seq = src->seq + start;
        if (src->q)
            target->q = src->q + start;
        target->len = len;
        target->flags = (target->flags | NOALLOC);
    }
}

size_t twobitseq(nucleotide* nucseqseq, size_t len, uint8_t** target) {
    size_t i;
    uint8_t* result;
    char* curbyte;
    size_t offset;
    size_t numel;
    
    numel = (len + 3) / 4;
    result = *target;
    result = realloc(result,numel*sizeof(uint8_t));
    /*note that if len is not a multiple of 4, the last byte will have A's appended to it*/
    curbyte = result;
    offset = 6;
    *curbyte = 0;
    for (i = 0;i < len;i++) {
        *curbyte |= (nucseqseq[i]<<offset);
        if (offset == 0) {
            offset = 6;
            curbyte++;
            *curbyte = 0;
        }
        else {
            offset -= 2;
        }
    }
    *target = result;
    return numel;
}
size_t twobitrcseq(nucleotide* nucseqseq, size_t len, uint8_t** target) {
    size_t i;
    uint8_t* result;
    char* curbyte;
    size_t offset;
    size_t numel;
    nucleotide n;

    numel = (len + 3) / 4;
    result = *target;
    result = realloc(result, numel * sizeof(uint8_t));
    /*note that if len is not a multiple of 4, the last byte will have A's appended to it*/
    curbyte = result;
    offset = 6;
    *curbyte = 0;
    for (i = 1;i <= len;i++) {
        n = nuc_complement(nucseqseq[len-i]);
        *curbyte |= (n << offset);
        if (offset == 0) {
            offset = 6;
            curbyte++;
            *curbyte = 0;
        }
        else {
            offset -= 2;
        }
    }
    *target = result;
    return numel;
}

nucseq** nucseq_winread(nucseq* src, size_t* count, size_t winsize, size_t minstep) {
    size_t seqlen;
    size_t nnewseqs;
    size_t i;
    size_t seqnamelen;
    size_t newseqnamelen;
    nucseq** result;
    double unistep;
    if (!src)return NULL;
    seqlen = src->len;
    if (winsize == 0)winsize = 1000;
    if (minstep == 0)minstep = winsize;
    if (seqlen < winsize) {
        *count = 0;
        return NULL;
    }
    if (src->name)seqnamelen = strlen(src->name);
    else seqnamelen = 0;
    if (*count == 0) {
        /*get as many sequences as possible*/
        nnewseqs = (seqlen - winsize) / minstep + 1;
        *count = nnewseqs;
        newseqnamelen = seqnamelen + (size_t)log10((double)nnewseqs) + 2;
        result = (nucseq**)malloc(sizeof(nucseq*)*(nnewseqs + 1));
        for (i = 0;i < nnewseqs;i++) {
            result[i] = (nucseq*)malloc(sizeof(nucseq));
            force_clear_nucseq(result[i]);
        }
        result[i] = NULL;
        for (i = 0;i < nnewseqs;i++) {
            subseq(result[i], src, i*minstep, winsize);
            if (seqnamelen) {
                result[i]->name = malloc(newseqnamelen + 1);
#ifdef _WIN32
                sprintf_s(result[i]->name, newseqnamelen + 1, "%s.%lld", src->name, i);
#else
                sprintf(result[i]->name, "%s.%zd", src->name, i);
#endif
            }
        }
    }
    else {
        /*it appears that the user has provided the amount of desired output sequences*/
        nnewseqs = (seqlen - winsize) / minstep + 1;
        if (*count > nnewseqs) *count = nnewseqs; /* if the user asks for too many sequences, change the amount*/
        else nnewseqs = *count;
        unistep = (double)(seqlen - winsize - 1) / (double)minstep;
        result = (nucseq**)malloc(sizeof(nucseq*)*(nnewseqs + 1));
        newseqnamelen = seqnamelen + (size_t)log10((double)nnewseqs) + 2;
        result = (nucseq**)malloc(sizeof(nucseq*)*(nnewseqs + 1));
        for (i = 0;i < nnewseqs;i++) {
            result[i] = (nucseq*)malloc(sizeof(nucseq));
            force_clear_nucseq(result[i]);
        }
        result[i] = NULL;
        /* Ideally, sequences should be equally distributed across the whole genome*/
        for (i = 0;i < nnewseqs;i++) {
            subseq(result[i], src, (size_t)(i*unistep + 0.01), winsize);
            if (seqnamelen) {
                result[i]->name = malloc(newseqnamelen + 1);
#ifdef _WIN32
                sprintf_s(result[i]->name, newseqnamelen + 1, "%s.%lld", src->name, i);
#else
                sprintf(result[i]->name, "%s.%zd", src->name, i);
#endif
            }
        }
    }
    return result;
}
nucseq** nucseq_winreadmultiple(nucseq** src, size_t nsequences, size_t* output, size_t winsize, size_t minstep) {
    nucseq** result;
    nucseq** tmpresults;
    size_t outsequences;
    size_t totnumseqs;
    size_t i;
    *output = 0;
    totnumseqs = 0;
    result = NULL;
    if (nsequences > 0) {
        for (i = 0;i < nsequences;i++) {
            if (src[i]) {
                outsequences = 0;
                tmpresults = nucseq_winread(src[i], &outsequences, winsize, minstep);
                if (tmpresults != NULL) {
                    totnumseqs += outsequences;
                    if (!result)
                        result = tmpresults;
                    else {
                        result = (nucseq**)realloc(result, sizeof(nucseq*)*(totnumseqs + 1));
                        memcpy(result + (totnumseqs - outsequences), tmpresults, sizeof(nucseq*)*(outsequences + 1));
                        free(tmpresults);
                    }
                }
            }
        }
    }
    else {
        i = 0;
        while (src[i]) {
            outsequences = 0;
            tmpresults = nucseq_winread(src[i], &outsequences, winsize, minstep);
            if (tmpresults != NULL) {
                totnumseqs += outsequences;
                if (!result)
                    result = tmpresults;
                else {
                    result = (nucseq**)realloc(result, sizeof(nucseq*)*(totnumseqs + 1));
                    memcpy(result + (totnumseqs - outsequences), tmpresults, sizeof(nucseq*)*(outsequences + 1));
                    free(tmpresults);
                }
            }
            i++;
        }
    }
    *output = totnumseqs;
    return result;
}

nucseq** nucseq_cutout(nucseq* src, nucseq* open, nucseq* close, size_t maxlen, size_t* output, uint32_t flags) {
    size_t begin, end, bfind, efind;
    int hasbegin, hasend;
    int reqbegin, reqend, allowoverlap, considerrc;
    int is_acceptable;
    int beginstatus, endstatus;
    nucseq** result;
    nucseq** tmpresult;
    nucseq tmprevopen, tmprevclose;
    size_t nalloc, newalloc;
    size_t nresult;
    size_t tmpn;
    size_t newlen;
    begin = end = bfind = efind = 0;
    hasbegin = hasend = 0;

    if (close && close->len <= 0) close = NULL;
    if (open && open->len <= 0) open = NULL;
    
    reqbegin = ((flags & NUCCUT_REQUIRE_BEGIN) != 0);
    reqend = ((flags & NUCCUT_REQUIRE_END) != 0);
    allowoverlap = ((flags & NUCCUT_ALLOWOVERLAP) != 0);
    considerrc = ((flags & NUCCUT_CONSIDERRC) != 0);

    bfind = nucsearch(src, open);
    if (!close && !reqend) {
        if (maxlen > 0)
            efind = bfind + maxlen;
        else
            efind = src->len;
    }
    else
        efind = nucsearch(src, close);

    nalloc = 16;
    nresult = 0;
    result = calloc(nalloc, sizeof(nucseq*));

    beginstatus = (open != NULL && open->len > 0 && bfind < src->len);
    endstatus = (close != NULL && close->len > 0 && efind < src->len);
    if (reqbegin && (!open || open->len < 1)) beginstatus = endstatus = 0;
    if (reqend && (!close || close->len < 1)) beginstatus = endstatus = 0;
    while ( beginstatus || endstatus) {
        if (bfind < src->len && efind > bfind) {
            begin = bfind;
            hasbegin = 1;
            if (allowoverlap)
                bfind = nucsearch_from(src, open, begin + 1);
            else
                bfind = nucsearch_from(src, open, efind + 1);
            if (bfind == 0)bfind = src->len;
        }
        if (efind < src->len) {
            end = efind;
            if (allowoverlap)
                efind = nucsearch_from(src, close, end + 1);
            else
                efind = nucsearch_from(src, close, bfind + 1);
            hasend = 1;
        }
        else {
            hasend = 0;
            end = src->len;
            if (close)end -= close->len;
        }
        is_acceptable = ((hasbegin || !reqbegin) && (hasend || !reqend));
        if (is_acceptable) {
            newlen = end - begin;
            if(close) newlen += close->len;
            if (newlen > maxlen && maxlen > 0) {
                result[nresult] = (nucseq*)calloc(1, sizeof(nucseq));
                if (!reqend && hasbegin)
                    subseq(result[nresult], src, begin, maxlen);
                else if (!reqbegin && hasend)
                    subseq(result[nresult], src, end + close->len - maxlen, maxlen);
                else {
                    is_acceptable = 0;
                    free(result[nresult]);
                }
            }
            else {
                result[nresult] = (nucseq*)calloc(1, sizeof(nucseq));
                subseq(result[nresult], src, begin, newlen);
            }
            if (is_acceptable) {
                nresult++;
                if (nresult >= nalloc) {
                    nalloc *= 2;
                    result = (nucseq**)realloc(result, sizeof(nucseq*)*nalloc);
                    for (tmpn = nresult;tmpn < nalloc;tmpn++) {
                        result[tmpn] = (nucseq*)calloc(1,sizeof(nucseq));
                    }
                }
            }
        }
        beginstatus = (open != NULL && open->len > 0 && bfind < src->len);
        endstatus = (close != NULL && close->len > 0 && efind < src->len);
    }
    if (considerrc) {
        flags = 0;
        force_clear_nucseq(&tmprevclose);
        force_clear_nucseq(&tmprevopen);
        if (reqbegin) flags |= NUCCUT_REQUIRE_END;
        if (reqend) flags |= NUCCUT_REQUIRE_BEGIN;
        if (allowoverlap) flags |= NUCCUT_ALLOWOVERLAP;
        nucseqrevcomp(open, &tmprevopen);
        nucseqrevcomp(close, &tmprevclose);
        tmpresult = nucseq_cutout(src, &tmprevclose, &tmprevopen, maxlen, &newalloc, flags);
        result = (nucseq**)realloc(result, sizeof(nucseq*)*(nresult+newalloc));
        for (tmpn = 0;tmpn < newalloc;tmpn++) {
            result[nresult + tmpn] = calloc(1, sizeof(nucseq));
            nucseqrevcomp(tmpresult[tmpn], result[nresult + tmpn]);
            clear_nucseq(tmpresult[tmpn]);
            free(tmpresult[tmpn]);
        }
        nresult += newalloc;
        clear_nucseq(&tmprevopen);
        clear_nucseq(&tmprevclose);
        free(tmpresult);
    }
    else
        result = (nucseq**)realloc(result, sizeof(nucseq*)*nresult);
    *output = nresult;
    return result;
}
nucseq** nucseq_cutoutmultiple(nucseq** src, size_t nsequences, nucseq* open, nucseq* close, size_t maxlen, size_t* output, uint32_t flags) {
    size_t osize,osizefull;
    size_t i,j;
    nucseq** result = NULL;
    nucseq** tmpresult;
    osizefull = osize = 0;
    for (i = 0;i < nsequences;i++) {
        tmpresult = nucseq_cutout(src[i], open, close, maxlen, &osize, flags);
        result = (nucseq**) realloc(result, sizeof(nucseq*)*(osizefull+osize));
        for (j = 0;j < osize;j++) {
            result[j + osizefull] = tmpresult[j];
        }
        osizefull += osize;
        free(tmpresult);
    }
    *output = osizefull;
    return result;
}

uint32_t* oligocount(nucseq* src, int k_len)
{
    uint32_t* counts = NULL;
    int64_t curn = 0;
    uint32_t numvals = 1;
    size_t contiglen = 0;
    int64_t i_;
    int64_t j_;
    int64_t k_;
    // high lengths do not make sense with this implementation
    if (k_len > 0 && k_len <= MAXOLIGOSIG) {
        // count the number of cells we need to have. No point in using pow here.
        for (curn = 0;curn < k_len;curn++) numvals *= 4;
        counts = (uint32_t*)malloc(sizeof(uint32_t)*numvals);
        if (!counts) {
            fprintf(stderr, "Not enough memory for count");
        }
        memset(counts, 0, sizeof(uint32_t)*numvals);
        for (curn = 0;curn < (int64_t)(src->len) - k_len + 1;curn++) {
            if (src->seq[curn] >= 0)
            {
                contiglen++;
                i_ = 0;
                j_ = k_len << 1;
                for (k_ = curn; k_ < curn + k_len; k_++)
                {
                    j_ -= 2;
                    i_ += src->seq[k_] << j_;
                }
                if (i_<numvals && i_ >= 0)
                    counts[i_]++;
            }
            else
            {
                contiglen = 0;
            }
        }
    }
    return counts;
}
uint32_t* oligocount_2strand(nucseq* src, int k_len)
{
    uint32_t* counts = NULL;
    int64_t curn = 0;
    uint32_t numvals = 1;
    size_t contiglen = 0;
    int64_t i_;
    int64_t j_;
    int64_t k_;
    // high lengths do not make sense with this implementation
    if (k_len > 0 && k_len <= MAXOLIGOSIG) {
        // count the number of cells we need to have. No point in using pow here.
        for (curn = 0;curn < k_len;curn++) numvals *= 4;
        counts = (uint32_t*)malloc(sizeof(uint32_t)*numvals);
        if (!counts) {
            fprintf(stderr, "Not enough memory for count");
        }
        memset(counts, 0, sizeof(uint32_t)*numvals);
        for (curn = 0;curn < (int64_t)(src->len) - k_len + 1;curn++) {
            if (src->seq[curn] >= 0)
            {
                contiglen++;
                i_ = 0;
                j_ = k_len << 1;
                for (k_ = curn; k_ < curn + k_len; k_++)
                {
                    j_ -= 2;
                    i_ += src->seq[k_] << j_;
                }
                if (i_<numvals && i_ >= 0)
                    counts[i_]++;
                i_ = 0;
                j_ = k_len << 1;
                for (k_ = curn + k_len - 1; k_ >= curn; k_--)
                {
                    j_ -= 2;
                    i_ += nuc_complement(src->seq[k_]) << j_;
                }
                if (i_<numvals && i_ >= 0)
                    counts[i_]++;
            }
            else
            {
                contiglen = 0;
            }
        }
    }
    return counts;
}
uint64_t* oligocount64_2strand(nucseq* src, int k_len)
{
    uint64_t* counts = NULL;
    int64_t curn = 0;
    uint64_t numvals = 1;
    size_t contiglen = 0;
    int64_t i_;
    int64_t j_;
    int64_t k_;
    // high lengths do not make sense with this implementation
    if (k_len > 0 && k_len <= MAXOLIGOSIG) {
        // count the number of cells we need to have. No point in using pow here.
        for (curn = 0;curn < k_len;curn++) numvals *= 4;
        counts = (uint64_t*)malloc(sizeof(uint64_t)*numvals);
        if (!counts) {
            fprintf(stderr, "Not enough memory for count");
        }
        memset(counts, 0, sizeof(uint64_t)*numvals);
        for (curn = 0;curn < (int64_t)(src->len) - k_len + 1;curn++) {
            if (src->seq[curn] >= 0)
            {
                contiglen++;
                i_ = 0;
                j_ = k_len << 1;
                for (k_ = curn; k_ < curn + k_len; k_++)
                {
                    j_ -= 2;
                    i_ += src->seq[k_] << j_;
                }
                if ((uint64_t)i_<numvals && i_ >= 0)
                    counts[i_]++;
                i_ = 0;
                j_ = k_len << 1;
                for (k_ = curn + k_len - 1; k_ >= curn; k_--)
                {
                    j_ -= 2;
                    i_ += nuc_complement(src->seq[k_]) << j_;
                }
                if ((uint64_t)i_<numvals && i_ >= 0)
                    counts[i_]++;
            }
            else
            {
                contiglen = 0;
            }
        }
    }
    return counts;
}

double nucleotide_identity(nucseq* A, nucseq* B) {
    size_t i;
    size_t maxi;
    size_t result_i;
    double result;
    maxi = (A->len > B->len ? B->len : A->len);
    result_i = 0;
    for (i = 0;i < maxi;i++) {
        if (A->seq[i] != B->seq[i])result_i++;
    }
    result = (double)(result_i) / (double)(maxi);
    return 1-result;
}
double amino_acid_identity(nucseq* A, nucseq* B) {
    size_t i;
    size_t maxi;
    size_t result_i;
    char A_aa;
    char B_aa;
    double result;
    maxi = (A->len > B->len ? B->len : A->len);
    maxi = (maxi - (maxi % 3));
    result_i = 0;
    for (i = 0;i < maxi;i+=3) {
        A_aa = translate_nuc(A->seq[i], A->seq[i + 1], A->seq[i + 2]);
        B_aa = translate_nuc(B->seq[i], B->seq[i + 1], B->seq[i + 2]);
        if (A_aa != B_aa)result_i++;
    }
    result = (double)(result_i) / (double)(maxi/3);
    return 1-result;
}
double nucseq_GC(nucseq * src)
{
    size_t numnuc = 0;
    size_t numGC = 0;
    size_t i_ = 0;
    for (i_ = 0;i_ < src->len;i_++)
    {
        if (src->seq[i_] == nucC || src->seq[i_] == nucG)
            numGC++;
        if (src->seq[i_] != nucN)
            numnuc++;
    }
    return ((double)numGC) / ((double)numnuc);
}
double* nucseq_GC_multiple(nucseq** src, size_t numseq)
{
    size_t i;
    double* result;
    result = (double*) malloc(numseq*sizeof(double));
    for (i = 0;i < numseq;i++) {
        result[i] = nucseq_GC(src[i]);
    }
    return result;
}
size_t minimal_qscore(nucseq * seq)
{
    int result = 0;
    size_t i = 0;
    if (seq && seq->q && seq->len>0) {
        result = seq->q[0];
        while (++i < seq->len)
            result = (result > seq->q[i]) ? seq->q[i] : result;
    }
    return result;
}

double* freqsig(uint32_t* counts, int k_len)
{
    double* result = NULL;
    size_t curn = 0;
    uint32_t numvals = 1;
    size_t totcount = 0;
    size_t klenl = (size_t)k_len;
    // high lengths do not make sense with this implementation
    if (k_len >= 0 && k_len <= MAXOLIGOSIG) {
        // count the number of cells we need to have. No point in using pow here.
        for (curn = 0;curn < klenl;curn++) numvals *= 4;
        result = (double*)malloc(sizeof(double)*numvals);
        if (!result) {
            fprintf(stderr, "Not enough memory for frequency storage");
        }
        for (curn = 0;curn < numvals;curn++)
            totcount += counts[curn];
        for (curn = 0;curn < numvals;curn++)
            result[curn] = ((double)counts[curn]) / ((double)totcount);
    }
    return result;
}

inline int count_set_bits(uint64_t value) {
    // For other ways to count bits, refer to:
    //  http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable 
    uint64_t v = value;
    unsigned int c;
    for (c = 0; v; c++) {
        v &= v - 1; // clear the least significant bit set
    }
    return c;
}

uint64_t bindex(uint64_t S, uint64_t B, int setbits) {
    // S represents [X_1 X_2 ... X_k] as SUM_n(X_n*4^n)
    // B represents [B_1 B_2 ... B_k] as SUM_n(B_n*2^n)
    // setbits is equal to SUM_n(B_n)
    int c, X;
    uint64_t mask = 1;
    uint64_t div = 1;
    uint64_t mult = 1;
    uint64_t toret = -1;
    if (B != 0) {
        toret = 0;
        for (c = 0;c < setbits;c++) {
            // count up to the next next non-zero bit
            while ((mask&B) == 0) {
                mask <<= 1;
                div *= 4;
            }
            // get the appropriate digit
            X = (S / div) % 4;
            toret += X*mult;
            // go to the next digit
            mult *= 4;
            div *= 4;
            mask <<= 1;
        }
    }
    return toret;
}

double* karlinsig(uint32_t* counts, int k)
{
    uint64_t idx, i;
    double* rho;
    double value;
    uint64_t s, s_max;
    uint64_t b, b_max;
    int* nCkSB;
    uint64_t** CkSB = NULL;
    uint64_t CkS0 = 0;
    int match = 0;
    int k_parity = k % 2;
    int denoms = 0;
    int zero_level = 0;
    // first calculate the number of elements in counts
    // number of S combinations: 4^k
    s_max = (uint64_t)pow(4, k);
    // number of B combinations: 2^k-2 ( [0 0 ... 0] is not valid and [1 1 ... 1] is already stored in counts)
    b_max = ((uint64_t)1 << k) - 2;
    // allocate rho
    rho = (double*)malloc(sizeof(double)*s_max);
    // allocate CkSB
    CkSB = (uint64_t**)malloc(sizeof(uint64_t)*b_max);
    nCkSB = (int*)malloc(sizeof(int)*b_max);
    for (b = 0;b < b_max;b++) {
        nCkSB[b] = count_set_bits(b + 1); // 4^(count_set_bits(b+1))
        i = 1ULL << (2 * nCkSB[b]);
        CkSB[b] = (uint64_t*)malloc(sizeof(uint64_t)*i);
        memset(CkSB[b], 0, sizeof(uint64_t)*i);
    }

    // Calculate CkSB
    // iterate over all elements
    for (s = 0; s < s_max; s++) {
        // for each combination add the value to the proper place
        for (b = 0;b < b_max;b++) {
            idx = bindex(s, b + 1, nCkSB[b]);
            CkSB[b][idx] += counts[s];
        }
        CkS0 += counts[s];
    }
    // calculate rho
    for (s = 0;s < s_max;s++) {
        if (counts[s] > 0)
            rho[s] = (double)counts[s] / (double)CkS0;
        else
            rho[s] = 1.0 / (double)CkS0 / (double)s_max;
        for (b = 0;b < b_max;b++) {
            idx = bindex(s, b + 1, nCkSB[b]);
            if (CkSB[b][idx] > 0)
                value = (double)CkSB[b][idx] / (double)CkS0;
            else
                value = 1 / (double)CkS0 / (double)nCkSB[b];
            if (nCkSB[b] % 2 == k_parity)
                rho[s] *= value;
            else
                rho[s] /= value;
        }
    }
    // free CkSB and nCkSB - we do not need them anymore
    if (CkSB) {
        for (b = 0;b < b_max;b++) {
            if (CkSB[b])
                free(CkSB[b]);
        }
        free(CkSB);
    }
    if (nCkSB)
        free(nCkSB);
    return rho;
}

void random_nucseq(nucseq* target, size_t nucseq_len, uint32_t seed) {
    size_t i;
    clear_nucseq(target);
    target->seq = malloc(sizeof(nucleotide)*nucseq_len);
    target->len = nucseq_len;
    if (seed < 0x10000)
        seed = seed + 0xffff;
    for (i = 0;i < nucseq_len;i++) {
        seed = rand_u32_minstd(seed);
        target->seq[i] = seed % 4;
    }
}
void random_weighed_nucseq(nucseq* target, size_t nucseq_len, uint32_t seed, double* frequencies, int kmer_len) {
    size_t i,j;
    size_t kindex;
    size_t veclen, pkp, len1; /*pkp = previous key possibilities*/
    uint32_t seed1;
    double tmp[4];
    double* probabilities;
    double tmpf;
    clear_nucseq(target);
    target->seq = malloc(sizeof(nucleotide)*nucseq_len);
    target->len = nucseq_len;
    veclen = 1LL << (2 * kmer_len);
    pkp = veclen >> 2;
    probabilities = malloc(sizeof(double)*veclen);
    if (seed < 0x10000)
        seed = seed + 0xffff;
    seed1 = seed;
    if (frequencies) {
        for (i = 0;i < pkp;i++) {
            tmpf = 0;
            for (j = (i << 2);j < (i << 2) + 4;j++)
                tmpf += frequencies[j];
            for (j = (i << 2);j < (i << 2) + 4;j++)
                probabilities[j] = (frequencies[j] / tmpf)*32767.0;
        }
    }
    else {
        for (i = 0;i < pkp;i++) {
            for (j = 0;j < 4;j++) {
                seed = rand_u32_minstd(seed);
                tmp[j] = (double)(seed % 32768);
            }
            tmpf = tmp[0] + tmp[1] + tmp[2] + tmp[3];
            for (j = 0;j < 4;j++)
                probabilities[j + (i << 2)] = (tmp[j] / tmpf) * 32767.0;
        }
    }
    seed = seed1;
    len1 = kmer_len-1;
    if (len1 > nucseq_len)len1 = nucseq_len;
    kindex = 0;
    for (i = 0;i < len1;i++) {
        seed = rand_u32_minstd(seed);
        target->seq[i] = seed % 4;
        kindex = (kindex | target->seq[i]);
        kindex = (kindex << 2);
    }
    for (i = len1;i < nucseq_len;i++) {
        seed = rand_u32_minstd(seed);
        tmpf = (double)(seed%32768);
        j = kindex;
        while (tmpf > probabilities[j]) {
            tmpf -= probabilities[j];
            j++;
        }
        target->seq[i] = (nucleotide)(j-kindex);
        kindex = (kindex | j) % pkp;
        kindex = (kindex << 2);
    }
}

void nucseq_oligocount_to_ht(nucseq* sequence, ht64_t* target, int k) {
    size_t i;
    size_t maxi;
    size_t stk;
    uint8_t* tbs;
    int nf;
    size_t ignorecount;
    size_t keylen;
    stk = (size_t)k;
    maxi = sequence->len - stk + 1;
    tbs = NULL;
    ignorecount = 0;
    for (i = 0;i < stk;i++) {
        if (sequence->seq[i] == nucN)
            ignorecount = i+1;
    }
    for (i = 0;i < maxi;i++) {
        if (sequence->seq[i+stk-1] == nucN) {
            ignorecount = k;
        }
        if (ignorecount > 0)
            ignorecount--;
        else {
            keylen = twobitseq(sequence->seq + i, stk, &tbs);
            ht64_inc(target, tbs, stk, 1, &nf);
        }
    }
    free(tbs);
}
ht64_t* kmer_table(nucseq** allsequences, size_t nseqs, int k) {
    /* calculate total size */
    uint64_t nbins;
    ht64_t* result;
    size_t i;
    nbins = 0;
    for (i = 0;i < nseqs;i++) {
        nbins += allsequences[i]->len;
    }
    if (nbins > 20000000)nbins = 20000000;
    result = ht64_alloc_size(nbins);
    for (i = 0;i < nseqs;i++) {
        nucseq_oligocount_to_ht(allsequences[i], result, k);
    }
    return result;
}

void nucseq_oligocount_to_hashht(nucseq* sequence, ht64_t* target, int k, hashdesc_t* hash) {
    size_t i;
    size_t maxi;
    size_t stk;
    uint64_t hkey;
    uint8_t* tbs;
    uint8_t* tbsrc;
    int nf;
    size_t ignorecount;
    size_t keylen;
    stk = (size_t)k;
    maxi = sequence->len - stk + 1;
    tbs = NULL;
    tbsrc = NULL;
    ignorecount = 0;
    for (i = 0;i < stk;i++) {
        if (sequence->seq[i] == nucN)
            ignorecount = i + 1;
    }
    for (i = 0;i < maxi;i++) {
        if (sequence->seq[i + stk - 1] == nucN) {
            ignorecount = k;
        }
        if (ignorecount > 0)
            ignorecount--;
        else {
            keylen = twobitseq(sequence->seq + i, stk, &tbs);
            twobitrcseq(sequence->seq + i, stk, &tbsrc);
            if (memcmp(tbs, tbsrc, keylen) < 0) {
                hkey = hash_index(tbs, keylen, hash);
                ht64_set(target, tbs, keylen, hkey, &nf);
            }
            else {
                hkey = hash_index(tbsrc, keylen, hash);
                ht64_set(target, tbsrc, keylen, hkey, &nf);
            }
        }
    }
    free(tbs);
    free(tbsrc);
}
ht64_t* kmerhash_table(nucseq** allsequences, size_t nseqs, int k, hashdesc_t* hash) {
    uint64_t nbins;
    ht64_t* result;
    size_t i;
    nbins = 0;
    for (i = 0;i < nseqs;i++) {
        nbins += allsequences[i]->len;
    }
    if (nbins > 20000000)nbins = 20000000;
    result = ht64_alloc_size(nbins);
    for (i = 0;i < nseqs;i++) {
        nucseq_oligocount_to_hashht(allsequences[i], result, k, hash);
    }
    return result;
}

void nucseq_oligocount_to_minhashLht(nucseq* sequence, ht64_t* target, int k, hashdesc_t* hash, size_t mod) {
    size_t i;
    size_t maxi;
    size_t stk;
    uint64_t hkey;
    uint8_t* tbs;
    uint8_t* tbsrc;
    int nf;
    size_t ignorecount;
    size_t keylen;
    stk = (size_t)k;
    maxi = sequence->len - stk + 1;
    tbs = NULL;
    tbsrc = NULL;
    ignorecount = 0;
    for (i = 0;i < stk;i++) {
        if (sequence->seq[i] == nucN)
            ignorecount = i + 1;
    }
    for (i = 0;i < maxi;i++) {
        if (sequence->seq[i + stk - 1] == nucN) {
            ignorecount = k;
        }
        if (ignorecount > 0)
            ignorecount--;
        else {
            keylen = twobitseq(sequence->seq + i, stk, &tbs);
            twobitrcseq(sequence->seq + i, stk, &tbsrc);
            if (memcmp(tbs, tbsrc, keylen) < 0) {
                hkey = hash_index(tbs, keylen, hash);
                if (hkey%mod == 0)
                    ht64_set(target, tbs, keylen, hkey, &nf);
            }
            else {
                hkey = hash_index(tbsrc, keylen, hash);
                if (hkey%mod == 0)
                    ht64_set(target, tbsrc, keylen, hkey, &nf);
            }
        }
    }
    free(tbs);
    free(tbsrc);
}
ht64_t* kmerhash_minhashLtable(nucseq** allsequences, size_t nseqs, int k, hashdesc_t* hash, size_t mod) {
    uint64_t nbins;
    ht64_t* result;
    size_t i;
    nbins = 0;
    for (i = 0;i < nseqs;i++) {
        nbins += allsequences[i]->len;
    }
    if (nbins > 20000000)nbins = 20000000;
    result = ht64_alloc_size(nbins);
    for (i = 0;i < nseqs;i++) {
        nucseq_oligocount_to_minhashLht(allsequences[i], result, k, hash, mod);
    }
    return result;
}

size_t nucseq_minimizers(nucseq* sequence, size_t** target, int k, size_t win) {
    dlink64_t* nextbestlist;
    dlink64_t* tmp_dlink;
    dlink64_t* tmp_dlink2;
    dlink64_t** positions;
    size_t* resultarray;
    size_t actualen;
    size_t nalloc;
    size_t i, j, nmins, oldbest, ignorecount;
    nalloc = 16;
    if (sequence->len < win) {
        *target = NULL;
        return 0;
    }
    actualen = win - k + 1;
    nextbestlist = NULL;
    positions = (dlink64_t**)calloc(actualen, sizeof(dlink64_t*));
    /* parse the first window */
    ignorecount = 0;
    for (i = 0;i < k;i++) {
        if (sequence->seq[i] == nucN) ignorecount = i+1;
    }
    j = ignorecount;
    ignorecount = 0;
    for (i = j;i < win - k;i++) {
        tmp_dlink = nextbestlist;
        /* the list is sorted, so the element needs to be inserted at the proper place
           note: a heap might be better than a linked list... but it is a bit harder to implement
           */
        tmp_dlink2 = NULL;
        if (sequence->seq[i + k - 1] == nucN) ignorecount = k;
        while (tmp_dlink && (ignorecount > 0 || memcmp(sequence->seq + i, sequence->seq + tmp_dlink->value, k) > 0)) {
            tmp_dlink2 = tmp_dlink;
            tmp_dlink = tmp_dlink->next;
        }
        if (ignorecount <= 0) {
            if (tmp_dlink)
                positions[i] = dlink_insert_before(tmp_dlink, i);
            else
                positions[i] = dlink_insert_after(tmp_dlink2, i);
            if (tmp_dlink2 == NULL)
                nextbestlist = positions[i];
        }
        else {
            positions[i] = NULL;
            ignorecount--;
        }
    }
    resultarray = (size_t*)realloc(*target, sizeof(size_t)*nalloc);
    resultarray[0] = nextbestlist->value;
    oldbest = nextbestlist->value;
    nmins = 1;
    for (i = actualen;i < sequence->len - k + 1;i++) {
        j = i%actualen;
        if (positions[j] == nextbestlist)
            nextbestlist = positions[j]->next;
        if(positions[j])
            dlink_remove(positions[j]);
        tmp_dlink = nextbestlist;
        tmp_dlink2 = NULL;
        if (sequence->seq[i] == nucN) ignorecount = k;
        while (tmp_dlink && (ignorecount > 0 || memcmp(sequence->seq + i, sequence->seq + tmp_dlink->value, k) > 0)) {
            tmp_dlink2 = tmp_dlink;
            tmp_dlink = tmp_dlink->next;
        }
        if (ignorecount <= 0) {
            if (tmp_dlink)
                positions[j] = dlink_insert_before(tmp_dlink, i);
            else if (tmp_dlink2)
                positions[j] = dlink_insert_after(tmp_dlink2, i);
            if (tmp_dlink2 == NULL)
                nextbestlist = positions[j];
            if (oldbest != nextbestlist->value) {
                oldbest = nextbestlist->value;
                if (nmins >= nalloc) {
                    nalloc *= 2;
                    resultarray = (size_t*)realloc(resultarray, sizeof(size_t)*nalloc);
                }
                resultarray[nmins] = nextbestlist->value;
                nmins++;
            }
        }
        else{
            positions[j] = NULL;
            ignorecount--;
        }
    }
    for (j = 0;j < actualen;j++) {
        if (positions[j])
            dlink_remove(positions[j]);
    }
    resultarray = (size_t*)realloc(resultarray, sizeof(size_t)*nmins);
    *target = resultarray;
    ___added_values___(0);
    free(positions);
    return nmins;
}
void nucseq_minicount_to_hashht(nucseq* sequence, ht64_t* target, int k, size_t win, hashdesc_t* hash) {
    size_t i, nmins, minimizer, keylen;
    uint64_t hkey;
    size_t* minimizers;
    uint8_t* tbs, *tbsrc;
    int nf;
    size_t stk;
    stk = k;
    minimizers = NULL;
    nmins = nucseq_minimizers(sequence, &minimizers, k, win);
    tbs = tbsrc = NULL;
    for (i = 0;i < nmins;i++) {
        minimizer = minimizers[i];
        keylen = twobitseq(sequence->seq + i, stk, &tbs);
        twobitrcseq(sequence->seq + i, stk, &tbsrc);
        if (memcmp(tbs, tbsrc, keylen) < 0) {
            hkey = hash_index(tbs, keylen, hash);
            ht64_set(target, tbs, keylen, hkey, &nf);
        }
        else {
            hkey = hash_index(tbsrc, keylen, hash);
            ht64_set(target, tbsrc, keylen, hkey, &nf);
        }
    }
    free(minimizers);
    free(tbs);
    free(tbsrc);
}
ht64_t* minimizerhash_table(nucseq** allsequences, size_t nseqs, int k, size_t win, hashdesc_t* hash){
    uint64_t nbins;
    ht64_t* result;
    size_t i;
    nbins = 0;
    for (i = 0;i < nseqs;i++) {
        nbins += allsequences[i]->len;
    }
    if (nbins > 2000000)nbins = 2000000;
    result = ht64_alloc_size(nbins);
    for (i = 0;i < nseqs;i++) {
        nucseq_minicount_to_hashht(allsequences[i], result, k, win, hash);
    }
    return result;
}
int64_t* minhash_Msig(nucseq** allsequences, size_t nseqs, int k, size_t siglen, size_t genomesize_est) {
    ht64_t* ht_seq;
    hashdesc_t* hash;
    int64_t* valuetable;
    int64_t* valuetable2;
    size_t nkmers;
    size_t i, j;
    size_t mod;
    hash = hashdesc_alloc();
    mod = (genomesize_est / siglen) / 3 * 2;
    hashdesc_init_fingerprint64(hash);
    ht_seq = kmerhash_minhashLtable(allsequences, nseqs, k, hash, mod);
    nkmers = ht64_astables(ht_seq, NULL, NULL, &valuetable);
    vec_sorti64(valuetable, nkmers);
    valuetable2 = (int64_t*)malloc(sizeof(int64_t)*(siglen + 1));
    j = 0;
    i = 0;
    valuetable2[0] = nkmers;
    while (i < nkmers && j < siglen) {
        valuetable2[j + 1] = valuetable[i];
        j++;
        i++;
    }
    free(valuetable);
    if (siglen > j) {
        for (i = j;i < siglen;i++)
            valuetable2[i] = 0x7FFFFFFFFFFFFFFF;
        /* by using this number, we guarantee that these values will have the lowest priority*/
    }
    /*cleanup*/
    ht64_free(ht_seq);
    hashdesc_free(hash);
    return valuetable2;
}


double* multifreqsig(nucseq** allsequences, size_t nseqs, int k, size_t winsize, char* primer) {
    nucseq** subseqs;
    uint32_t* oligos;
    size_t maxvecid;
    size_t numel_sig;
    size_t numel_tot;
    size_t i, j, m, bestseqid;
    size_t nsubseqs;
    double** freqsigs;
    double* result;
    double* currentbestsig;
    double score;
    double bestscore;
    nucseq tmpseq = EMPTYSEQ;
    nucseq targetseq = EMPTYSEQ;
    nucseq_from_string(&tmpseq, primer);
    numel_sig = 1LL << (2 * k);
    numel_tot = numel_sig * numel_sig;
    maxvecid = 0;
    result = (double*)malloc(sizeof(double)*numel_tot);
    subseqs = nucseq_cutoutmultiple(allsequences, nseqs, &tmpseq, NULL, winsize-1+tmpseq.len, &nsubseqs, NUCCUT_ALLOWOVERLAP|NUCCUT_REQUIRE_BEGIN|NUCCUT_CONSIDERRC);
    if (nsubseqs == 0) {
        memset(result, 0, sizeof(double)*numel_tot);
        return result; /* ABORT */
    }
    freqsigs = (double**)malloc(sizeof(double*)*nsubseqs);
    for (j = 0;j < nsubseqs;j++) {
        oligos = oligocount_2strand(subseqs[j], k);
        freqsigs[j] = freqsig(oligos, k);
        free(oligos);
    }
    for (i = 0;i < numel_sig;i++) {
        random_weighed_nucseq(&targetseq, 3000, (uint32_t)i, NULL, 4);
        bestseqid = -1;
        bestscore = 0;
        currentbestsig = freqsigs[0];
        for (j = 0;j < nsubseqs;j++) {
            if (subseqs[j]->len == winsize-1+tmpseq.len) {
                score = amino_acid_identity(&targetseq,subseqs[j]);
                if (score > bestscore) {
                    currentbestsig = freqsigs[j];
                    bestscore = score;
                    bestseqid = j;
                }
            }
        }
        m = i*numel_sig;
        for (j = 0;j < numel_sig;j++) {
            result[j + m] = currentbestsig[j];
        }
    }
    clear_nucseq(&targetseq);
    for (j = 0;j < nsubseqs;j++) {
        free(freqsigs[j]);
    }
    free(freqsigs);
    nucseq_array_clear(subseqs, nsubseqs);
    clear_nucseq(&tmpseq);
    return result;
}
double* multikarlsig(nucseq** allsequences, size_t nseqs, int k, size_t winsize, size_t minstep) {
    nucseq** subseqs;
    uint32_t* oligos;
    size_t maxvecid;
    size_t outputamnt;
    size_t numel_sig;
    size_t numel_tot;
    size_t i, j, m;
    double* result;
    double* currentbestsig;
    double* tmpsig;
    numel_sig = 1LL << (2 * k);
    numel_tot = numel_sig * numel_sig;
    subseqs = nucseq_winreadmultiple(allsequences, nseqs, &outputamnt, winsize, minstep);
    maxvecid = 0;
    result = (double*)malloc(sizeof(double)*numel_tot);
    for (i = 0;i < numel_sig;i++) {
        oligos = oligocount(subseqs[0], k);
        currentbestsig = karlinsig(oligos, k);
        free(oligos);
        for (j = 1;j < outputamnt;j++) {
            oligos = oligocount(subseqs[0], k);
            tmpsig = karlinsig(oligos, k);
            free(oligos);
            if (tmpsig[i] > currentbestsig[i]) {
                free(currentbestsig);
                currentbestsig = tmpsig;
            }
            else
                free(tmpsig);
        }
        m = i*numel_sig;
        for (j = 0;j < numel_sig;j++) {
            result[j + m] = currentbestsig[j];
        }
        free(currentbestsig);
    }
    return result;
}

double* fullmarkovsig(uint32_t* counts, int k) {
    uint64_t full_max, leftright_max, centre_max;
    uint64_t i;
    uint64_t fullcount_i;
    double* result;
    uint64_t* left_i;
    uint64_t* right_i;
    uint64_t* centre_i;
    uint64_t* left_target;
    uint64_t* right_target;
    uint64_t* centre_target;
    double fullcount;
    double expectation;
    double tmpleft;
    double tmpright;
    double tmpcentre;
    uint64_t lmask;
    uint64_t rmask;
    uint64_t cmask;
    full_max = (uint64_t)pow(4, k);
    result = (double*)calloc(full_max, sizeof(double));
    /* if the intersection between left and right would be "shorter than 1", abort*/
    if (k < 3) return result;
    leftright_max = (full_max >> 2);
    centre_max = (full_max >> 4);
    
    lmask = ((full_max - 1) - 0b11);
    rmask = ((full_max - 1) >> 2);
    cmask = (lmask & rmask);

    left_i = (uint64_t*)calloc(leftright_max, sizeof(uint64_t));
    right_i = (uint64_t*)calloc(leftright_max, sizeof(uint64_t));
    centre_i = (uint64_t*)calloc(centre_max, sizeof(uint64_t));

    left_target = (uint64_t*)malloc(full_max*sizeof(uint64_t));
    right_target = (uint64_t*)malloc(full_max*sizeof(uint64_t));
    centre_target = (uint64_t*)malloc(full_max*sizeof(uint64_t));

    fullcount_i = 0;
    for (i = 0;i < full_max;i++) {
        fullcount_i += counts[i];
        left_target[i] = ((i&lmask) >> 2);
        right_target[i] = (i&rmask);
        centre_target[i] = ((i&cmask) >> 2);
        left_i[left_target[i]] += counts[i];
        right_i[right_target[i]] += counts[i];
        centre_i[centre_target[i]] += counts[i];
    }
    fullcount = (double)fullcount_i;

    for (i = 0;i < full_max;i++) {
        tmpleft = (double)(left_i[left_target[i]]);
        tmpright = (double)(right_i[right_target[i]]);
        tmpcentre = (double)centre_i[left_target[i]];
        expectation = tmpleft*tmpright / tmpcentre;
        result[i] = (double)counts[i] / expectation;
    }
    return result;
}
double* fullmarkovsig_zscore(uint32_t* counts, int k) {
    uint64_t full_max, leftright_max, centre_max;
    uint64_t i;
    uint64_t fullcount_i;
    double* result;
    uint64_t* left_i;
    uint64_t* right_i;
    uint64_t* centre_i;
    uint64_t* left_target;
    uint64_t* right_target;
    uint64_t* centre_target;
    double fullcount;
    double expectation;
    double variance;
    double tmpleft;
    double tmpright;
    double tmpcentre;
    uint64_t lmask;
    uint64_t rmask;
    uint64_t cmask;
    full_max = (uint64_t)pow(4, k);
    result = (double*)calloc(full_max, sizeof(double));
    /* if the intersection between left and right would be "shorter than 1", abort*/
    if (k < 3) return result;
    leftright_max = (full_max >> 2);
    centre_max = (full_max >> 4);

    lmask = ((full_max - 1) - 0b11);
    rmask = ((full_max - 1) >> 2);
    cmask = (lmask & rmask);

    left_i = (uint64_t*)calloc(leftright_max, sizeof(uint64_t));
    right_i = (uint64_t*)calloc(leftright_max, sizeof(uint64_t));
    centre_i = (uint64_t*)calloc(centre_max, sizeof(uint64_t));

    left_target = (uint64_t*)malloc(full_max * sizeof(uint64_t));
    right_target = (uint64_t*)malloc(full_max * sizeof(uint64_t));
    centre_target = (uint64_t*)malloc(full_max * sizeof(uint64_t));

    fullcount_i = 0;
    for (i = 0;i < full_max;i++) {
        fullcount_i += counts[i];
        left_target[i] = ((i&lmask) >> 2);
        right_target[i] = (i&rmask);
        centre_target[i] = ((i&cmask) >> 2);
        left_i[left_target[i]] += counts[i];
        right_i[right_target[i]] += counts[i];
        centre_i[centre_target[i]] += counts[i];
    }
    fullcount = (double)fullcount_i;

    for (i = 0;i < full_max;i++) {
        tmpleft = (double)(left_i[left_target[i]]);
        tmpright = (double)(right_i[right_target[i]]);
        tmpcentre = (double)centre_i[left_target[i]];
        expectation = tmpleft*tmpright / tmpcentre;
        variance = expectation*(tmpcentre - tmpleft)*(tmpcentre - tmpright) / (tmpcentre*tmpcentre);
        result[i] = ((double)counts[i] - expectation)/(sqrt(variance));
    }
    return result;
}

double* TETRAsig(uint32_t* counts)
{
    int k = 4;
    uint64_t idx, i;
    double* rho;
    double E, var;
    uint64_t s, s_max;
    uint64_t b, b_max;
    int* nCkSB;
    uint64_t** CkSB = NULL;
    uint64_t CkS0 = 0;
    double N123, N234, N23;
    int match = 0;
    int k_parity = k % 2;
    int denoms = 0;
    int zero_level = 0;
    // first calculate the number of elements in counts
    // number of S combinations: 4^k
    s_max = (uint64_t)pow(4, k);
    // number of B combinations: 2^k-2 ( [0 0 ... 0] is not valid and [1 1 ... 1] is already stored in counts)
    b_max = ((uint64_t)1 << k) - 2;
    // allocate rho
    rho = (double*)malloc(sizeof(double)*s_max);
    // allocate CkSB
    CkSB = (uint64_t**)malloc(sizeof(uint64_t)*b_max);
    nCkSB = (int*)malloc(sizeof(int)*b_max);
    for (b = 0;b < b_max;b++) {
        nCkSB[b] = count_set_bits(b + 1); // 4^(count_set_bits(b+1))
        i = 1ULL << (2 * nCkSB[b]);
        CkSB[b] = (uint64_t*)malloc(sizeof(uint64_t)*i);
        memset(CkSB[b], 0, sizeof(uint64_t)*i);
    }

    // Calculate CkSB
    // iterate over all elements
    for (s = 0; s < s_max; s++) {
        // for each combination add the value to the proper place
        for (b = 0;b < b_max;b++) {
            idx = bindex(s, b + 1, nCkSB[b]);
            CkSB[b][idx] += counts[s];
        }
        CkS0 += counts[s];
    }
    // calculate rho
    for (s = 0;s < s_max;s++) {
        if (counts[s] > 0)
            rho[s] = (double)counts[s];
        else
            rho[s] = 1.0 / (double)CkS0;
        /* the only difference between this and karlin signatures */
        N123 = (double)CkSB[0b1110 - 1][bindex(s, 0b1110, nCkSB[0b1110 - 1])];
        N234 = (double)CkSB[0b0111 - 1][bindex(s, 0b0111, nCkSB[0b0111 - 1])];
        N23 = (double)CkSB[0b0110 - 1][bindex(s, 0b0110, nCkSB[0b0110 - 1])];
        E = N123*N234 / N23;
        var = E*(N23 - N123)*(N23 - N234) / (N23*N23);

        rho[s] = (rho[s] - E) / (sqrt(var));
    }
    // free CkSB and nCkSB - we do not need them anymore
    if (CkSB) {
        for (b = 0;b < b_max;b++) {
            if (CkSB[b])
                free(CkSB[b]);
        }
        free(CkSB);
    }
    if (nCkSB)
        free(nCkSB);
    return rho;
}
size_t alignmatch(nucseq* toalign, nucseq* src, int64_t offset)
{
    size_t i_ = 0;
    size_t j_ = 0;
    size_t match = 0;
    if (offset < 0)
        j_ = (size_t)(-offset);
    else
        i_ = (size_t)offset;
    while (i_ < src->len && j_ < toalign->len) {
        if (toalign->seq[j_] == src->seq[i_])
            match++;
        i_++;
        j_++;
    }
    return match;
}

int64_t best_perfect_align(nucseq* toalign, nucseq* src)
{
    size_t offset;
    size_t len;
    size_t minlen;
    size_t j_, k_;
    nucseq tmpseq = EMPTYSEQ;
    nucseq remseq = EMPTYSEQ;
    size_t match = 0;
    size_t bestmatch = 0;
    int64_t best_offset;
    size_t intern_offset;
    int found = 0;
    len = toalign->len / 2;
    best_offset = src->len;
    // the minimal search length is defined by how many copies of a subsequence are likely to
    // exist in the source. This is a function of src->len, but we use the value below for now
    // TODO: THIS IS A TEMPORARY MEASURE - once the function is know it needs to go here
    minlen = toalign->len / 10;
    if (minlen > 20) minlen = 20;
    while (!found && len > minlen)
    {
        intern_offset = 0;
        // check all substrings
        while (intern_offset + len <= toalign->len && !found) {
            subseq_nocpy(&tmpseq, toalign, intern_offset, len);
            offset = nucsearch(src, &tmpseq);
            found = (offset < src->len);
            intern_offset += len;
        }
        intern_offset -= len;
        if (found) {
            // we've found something, now we need to check if this place is the best match
            best_offset = (int64_t)offset - (int64_t)intern_offset;
            match = len;
            while (offset < src->len) {
                // extend backward
                j_ = intern_offset - 1;
                k_ = offset - 1;
                while (j_ >= 0 && k_ >= 0 && toalign->seq[j_] == src->seq[k_]) {
                    j_--;
                    k_--;
                }
                match += intern_offset - 1 - j_;
                // extend forward
                j_ = intern_offset + len;
                k_ = offset + len;
                while (j_ <toalign->len && k_<src->len && toalign->seq[j_] == src->seq[k_]) {
                    j_++;
                    k_++;
                }
                match += j_ - (intern_offset + len);
                if (match > bestmatch) {
                    bestmatch = match;
                    best_offset = (int64_t)offset - (int64_t)intern_offset;
                }
                // search for the next occurence
                offset++;
                subseq_nocpy(&remseq, src, offset, src->len - offset);
                offset = nucsearch(&remseq, &tmpseq) + offset;
            }
        }
        else
            len = len / 2; // we need to divide further
    }
    return best_offset;
}

size_t extseq(nucseq* target, nucseq* extension, uint32_t extension_flags, char* ext_args)
{
    int64_t target_offset;
    int64_t toa, tob;
    size_t newlen;
    char* consensus;
    char* newq;
    char targq, targn, extq, extn;
    nucseq prefix = EMPTYSEQ;
    nucseq suffix = EMPTYSEQ;
    newlen = 0;
    size_t tzero, ezero, i_;
    size_t suffixstart;
    int edgesonly;

    edgesonly = extension_flags&EDGESONLY;
    extension_flags &= (~EDGESONLY); // remove the EDGESONLY flag
    if (!extension_flags) {
        // find the place which has the longest perfect match
        if (edgesonly && target->len>2 * extension->len)
        {
            suffixstart = target->len - extension->len;
            subseq_nocpy(&prefix, target, 0, extension->len);
            subseq_nocpy(&suffix, target, suffixstart, extension->len);
            toa = best_perfect_align(extension, &prefix);
            tob = best_perfect_align(extension, &suffix);
            // Prefer the end part when alignments are the same.
            // This ensures that if both fail to find something, target_offset will
            // have the proper value
            if (alignmatch(extension, &prefix, toa) > alignmatch(extension, &suffix, tob))
                target_offset = toa;
            else
                target_offset = tob + suffixstart;
        }
        else
            target_offset = best_perfect_align(extension, target);
        if (target_offset == target->len) return target->len;
        if (target_offset > 0) {
            ezero = (size_t)target_offset;
            tzero = 0;
            newlen = extension->len + ezero;
            if (newlen<target->len) newlen = target->len;
        }
        else {
            ezero = 0;
            tzero = (size_t)(-target_offset);
            newlen = target->len + tzero;
            if (newlen<extension->len) newlen = extension->len;
        }
        consensus = (char*)malloc(sizeof(char)*newlen);
        newq = (char*)malloc(sizeof(char)*newlen);
        if (!consensus || !newq) {
            fprintf(stderr, "Not enough memory for extension");
        }
        memcpy(consensus + tzero, target->seq, target->len);
        if (target->q)
            memcpy(newq + tzero, target->q, target->len);
        else
            memset(newq + tzero, FASTA_Q, sizeof(char));
        // no weighing is applied, the highest q-score remains (no value = FASTA_Q)
        for (i_ = 0; i_ < newlen; i_++) {
            // if outside the overlap region, just copy values
            if (i_ >= ezero && i_ < extension->len + ezero &&
                i_ >= tzero && i_ < target->len + tzero)
            {
                extq = (extension->q) ? extension->q[i_ - ezero] : FASTA_Q;
                extn = extension->seq[i_ - ezero];
                targq = (target->q) ? target->q[i_ - tzero] : FASTA_Q;
                targn = target->seq[i_ - tzero];
                if (targn == extn) {
                    consensus[i_] = targn;
                    newq[i_] = targq + extq;
                    if (newq[i_] > 96) newq[i_] = 96;
                }
                else {
                    if (targq > extq) {
                        consensus[i_] = targn;
                        newq[i_] = targq - extq;
                    }
                    else {
                        consensus[i_] = extn;
                        newq[i_] = extq - targq;
                    }
                }
            }
            else if (i_ >= ezero && i_ < extension->len + ezero) {
                consensus[i_] = extension->seq[i_ - ezero];
                newq[i_] = (extension->q) ? extension->q[i_ - ezero] : FASTA_Q;
            }
        }
        clear_nucseq(target);
        target->seq = consensus;
        target->q = newq;
        target->len = newlen;
    }
    else {
        fprintf(stderr, "Extension Flags are not supported in this version");
    }
    return newlen;
}

size_t concactseqs(nucseq* target, nucseq* left, nucseq* right) {
    size_t newlen = 0;
    clear_nucseq(target);
    if (left && right) {
        target->flags = 0;
        target->len = newlen;
        newlen = left->len + right->len;
        if (left->q && right->q) {
            target->q = malloc(sizeof(char)*newlen);
            memcpy(target->q, left->q, left->len);
            memcpy(target->q + left->len, right->q, right->len);
        }
        target->seq = malloc(sizeof(char)*newlen);
        memcpy(target->seq, left->seq, left->len);
        memcpy(target->seq + left->len, right->seq, right->len);
    }
    else {
        if (left) {
            subseq(target, left, 0, left->len);
        }
        if (right) {
            subseq(target, right, 0, right->len);
        }
    }
    return newlen;
}
size_t appendtoseq(nucseq* target, nucseq* right) {
    size_t newlen = 0;
    size_t oldlen = 0;
    if (target && right) {
        if (target->len > 0) {
            oldlen = target->len;
            newlen = oldlen + right->len;
            target->len = newlen;
            if (target->q && right->q) {
                target->q = realloc(target->q, sizeof(char)*newlen);
                memcpy(target->q + oldlen, right->q, right->len);
            }
            target->seq = realloc(target->seq, sizeof(char)*newlen);
            memcpy(target->seq + oldlen, right->seq, right->len);
        }
        else
            subseq(target, right, 0, right->len);
    }

    return newlen;
}

/* ******************
   Boyer-Moore Search
   ******************
   adapted from Wikipedia */

void make_delta1(size_t *delta1, nucseq* pat) {
    size_t i;
    char* patseq;
    patseq = pat->seq;
    for (i = 0; i < 4; i++) {
        delta1[i] = pat->len;
    }
    for (i = 0; i < pat->len - 1; i++) {
        if (pat->seq[i] >= 0 && pat->seq[i]<4)
            delta1[pat->seq[i]] = pat->len - 1 - i;
    }
}

int is_prefix(nucseq* pat, size_t pos) {
    size_t i;
    size_t suffixlen;
    char* patseq;
    suffixlen = pat->len - pos;;
    patseq = pat->seq;
    /* could also use the strncmp() library function here */
    for (i = 0; i < suffixlen; i++) {
        if (patseq[i] != patseq[pos + i]) {
            return 0;
        }
    }
    return 1;
}

size_t suffix_length(nucseq* pat, size_t pos) {
    size_t i;
    for (i = 0; (pat->seq[pos - i] == pat->seq[pat->len - 1 - i]) && (i < pos); i++);
    return i;
}

void make_delta2(size_t* delta2, nucseq* pat) {
    long long p;
    long long patlen;
    size_t last_prefix_index;
    size_t slen;
    patlen = pat->len;
    last_prefix_index = patlen - 1;

    // first loop
    for (p = patlen - 1; p >= 0; p--) {
        if (is_prefix(pat, p + 1)) {
            last_prefix_index = p + 1;
        }
        delta2[p] = last_prefix_index + (patlen - 1 - p);
    }

    // second loop
    for (p = 0; p < patlen - 1; p++) {
        slen = suffix_length(pat, p);
        if (pat->seq[p - slen] != pat->seq[patlen - 1 - slen]) {
            delta2[patlen - 1 - slen] = patlen - 1 - p + slen;
        }
    }
}

size_t nucsearch_from(nucseq* data, nucseq* pattern, size_t start) {
    long long i, j;
    size_t delta1[4];
    size_t *delta2 = NULL;
    size_t patlen;
    char* pat;
    char* string;
    if (!data) return 0;
    if (!pattern) return data->len;
    patlen = pattern->len;
    pat = pattern->seq;
    string = data->seq;
    delta2 = malloc((patlen + 1) * sizeof(size_t));
    if (!delta2) {
        fprintf(stderr, "Not enough memory for nucsearch");
    }
    // The empty pattern must be considered specially
    if (patlen == 0) {
        free(delta2);
        return 0;
    }
    make_delta1(delta1, pattern);
    make_delta2(delta2, pattern);

    i = patlen - 1 + start;
    while (i>=0 && i < (long long)data->len) {
        j = patlen - 1;
        // using nucmatch instead of == ensures that we never end up on N
        while (j >= 0 && _nucmatch(string[i], pat[j])) {
            --i;
            --j;
        }
        if (j < 0) {
            free(delta2);
            return i + 1;
        }
        i += (delta1[string[i]] > delta2[j]) ? delta1[string[i]] : delta2[j];
    }
    free(delta2);
    return data->len;
}
size_t nucsearch(nucseq* data, nucseq* pattern) {
    return nucsearch_from(data, pattern, 0);
}