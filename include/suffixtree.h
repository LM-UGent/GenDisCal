#ifndef DEF_SUFFIXTREE
#define DEF_SUFFIXTREE

#include <stdint.h>
#include "osportstd.h"

typedef struct suffixtree_t suftree_t;

#define SUFTREECOPYTEXT     0b10000
#define SUFTREENOCOPY       0b100000
#define SUFTREE_DEFAULT     0
#define SUFTREE_NAIVE       0b00001
suftree_t* suftree_from(char* data, int option, size_t datalen);
void suftree_free(suftree_t* target);
void suftree_printf(suftree_t* target);
void suftree_printnewick(suftree_t* target);
void suftree_printbacklinks(suftree_t* tree);
char* suftree_getstrptr(suftree_t* tree, size_t* output_len);
void suftree_save(suftree_t* tree, PF_t* target);
suftree_t* suftree_load(PF_t* target);

suftree_t* suftree_intersection(suftree_t* A, suftree_t* B);
size_t suftree_search(suftree_t* target,char* data, size_t datalen);
char* suftree_firstsharedregions(suftree_t* reference, suftree_t* query, size_t minlen);
#define DATAMARK_BOTTOMSTARS 0
#define DATAMARK_INTERVALS  1
void suftree_printmarkeddata(char* marks, char* data, size_t datalen, suftree_t* tree, int outputformat);
char* suftree_approximatesearch(char* data, size_t datalen, suftree_t* tree, size_t seedalign, size_t minlen, size_t maxgaperror, int add_unaligned_lengths);
double suftree_roughalign(char* data, size_t datalen, suftree_t* tree, size_t minmatch, size_t smallgapsize, int trim_mismatches, size_t* alnlen);
size_t suftree_datalen(suftree_t* tree);


#define KEEP_MARKED 0
#define KEEP_UNMARKED 1
#define KEEP_UNMARKED_INTERVALS 2
char** marks_split(char* data, char* marks, size_t datalen, size_t** output_lens, size_t* output_amount, int tokeep);

#endif

