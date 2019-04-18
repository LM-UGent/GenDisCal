#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nucseq.h"
#include "vecops.h"

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
            result[curpos] = nuc2char(target->seq[i_]);
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

// ******************
// Boyer-Moore Search
// ******************
// adapted from Wikipedia

void make_delta1(size_t *delta1, nucseq* pat) {
    size_t i;
    char* patseq;
    patseq = pat->seq;
    for (i = 0; i < 4; i++) {
        delta1[i] = -1;
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
    // could also use the strncmp() library function here
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

size_t nucsearch(nucseq* data, nucseq* pattern) {
    long long i, j;
    size_t delta1[4];
    size_t *delta2 = NULL;
    size_t patlen;
    char* pat;
    char* string;
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

    i = patlen - 1;
    while (i < (long long)data->len) {
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