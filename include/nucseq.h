#ifndef DEF_NUCSEQ
#define DEF_NUCSEQ

#include "osportstd.h"


#define FASTA_Q     30

#define SINGLE_STRAND_SEQ   0b00001
#define DOUBLE_STRAND_SEQ   0b00000
#define CIRCULAR_SEQ    0b00010
#define LINEAR_SEQ      0b00000
#define NOALLOC         0b00100
#define READONLY        0b01000

typedef int8_t nucleotide;

#define nucA    0
#define nucC    1
#define nucG    2
#define nucT    3
#define nucN    -1

#define is_purine(nuc) (x==nucA || x==nucG)
#define is_pyrimidine(nuc) (x==nucC || x==nucT)
nucleotide nuc_complement(nucleotide nuc);
nucleotide char2nuc(char c);
char nuc2char(nucleotide c);


typedef struct nucseq {
    size_t len;
    nucleotide* seq;
    char* q;
    uint32_t flags;
    char* name;
}nucseq;
#define EMPTYSEQ    {0,NULL,NULL,0,NULL}

void clear_nucseq(nucseq* target);
size_t nucseq_from_string(nucseq* target, char* str);
char* nucseq_tritext(nucseq* target);
nucseq** nucseq_array_from_fasta(PF_t* f, size_t* OUT_count, int saveseqnames, size_t minlen);
char* nucseq2fasta(nucseq* target, const char* name);
void nucseq2fastafile(nucseq* target, PF_t* f, const char* name);
void nucseqrevcomp(nucseq* src, nucseq* target);

void subseq(nucseq* target, nucseq* src, size_t start, size_t len);
void subseq_nocpy(nucseq* target, nucseq* src, size_t start, size_t len);

nucseq** nucseq_winread(nucseq* src, size_t* count, size_t winsize, size_t minstep);
nucseq** nucseq_winreadmultiple(nucseq** src, size_t nsequences, size_t* output, size_t winsize, size_t minstep);
#define NUCCUT_REQUIRE_BEGIN    1
#define NUCCUT_REQUIRE_END      2
#define NUCCUT_ALLOWOVERLAP     4
#define NUCCUT_CONSIDERRC       8
nucseq** nucseq_cutout(nucseq* src, nucseq* open, nucseq* close, size_t maxlen, size_t* output, uint32_t flags);
nucseq** nucseq_cutoutmultiple(nucseq** src, size_t nsequences, nucseq* open, nucseq* close, size_t maxlen, size_t* output, uint32_t flags);
#define MAXOLIGOSIG     12
uint32_t* oligocount(nucseq* src, int k_len);
uint32_t* oligocount_2strand(nucseq* src, int k_len);
uint64_t* oligocount64_2strand(nucseq* src, int k_len);

double nucseq_GC(nucseq* src);
double* nucseq_GC_multiple(nucseq** src, size_t numseq);
size_t minimal_qscore(nucseq* seq);

double* freqsig(uint32_t* counts, int k_len);
double* TETRAsig(uint32_t* counts);
double* karlinsig(uint32_t* counts, int k_len);
int64_t* minhash_Msig(nucseq** allsequences, size_t nseqs, int k, size_t siglen, size_t genomesize_est);
double* multifreqsig(nucseq** allsequences, size_t nseqs, int k, size_t winsize, char* primer);
double* multikarlsig(nucseq** allsequences, size_t nseqs, int k, size_t winsize, size_t minstep);
double* fullmarkovsig(uint32_t* counts, int k);
double* fullmarkovsig_zscore(uint32_t* counts, int k);

size_t alignmatch(nucseq* toalign, nucseq* src, int64_t offset);
int64_t best_perfect_align(nucseq* target, nucseq* src);

#define PERFECT_MATCH   0
#define MISMATCH       -1
#define MATCH           1

#define ALLOW_GAPS      1
#define ALLOW_SNPS      2
#define EDGESONLY       4
size_t extseq(nucseq* target, nucseq* src, uint32_t extension_flags, char* ext_args);
size_t concactseqs(nucseq* target, nucseq* left, nucseq* right);
size_t appendtoseq(nucseq* target, nucseq* right);
size_t nucsearch(nucseq* data, nucseq* pattern);
size_t nucsearch_from(nucseq* data, nucseq* pattern, size_t start);
#endif
