#ifndef DEF_TAXONOMY
#define DEF_TAXONOMY

#include "osportstd.h"

typedef struct lineage_t lineage_t;
typedef struct taxtree_t taxtree_t;

void lineage_free(lineage_t* lineage);
int lineage_majornamerelatedness(lineage_t* A, lineage_t* B);
int64_t lineage_idbytaxrank(lineage_t* lineage, double level);
double lineage_parentlevel(lineage_t* lineage, double level);
int lineage_includes_id(lineage_t* lineage, int64_t id);


#define TAXFMT_PCOF_S   0
#define TAXFMT_GTDB     1
#define TAXFMT_NCBI     2
#define TAXFMT_CSV      0
#define TAXFMT_TSV      1

char** PCOF_S_majornames(char* str, size_t* outlen);

taxtree_t* taxtree_alloc();
taxtree_t* taxtree_fromlineages(PF_t* file, int fmt);
taxtree_t* taxtree_fromrelations(PF_t* file, int fmt);
void taxtree_addleaves(taxtree_t* tree, PF_t* file, int fmt);
void taxtree_rename(taxtree_t* tree, PF_t* file, int fmt);
void taxtree_free(taxtree_t* target);

int64_t taxtree_getid(taxtree_t* tree, char* name);
int64_t taxtree_getparentid(taxtree_t* tree, int64_t id);
lineage_t* taxtree_lineageof(taxtree_t* tree, int64_t id);
char* taxtree_getnameptr(taxtree_t* tree, int64_t id);
char* taxtree_getlabelptr(taxtree_t* tree, int64_t id);
char* taxtree_generatereadablename(taxtree_t* tree, int64_t id);
int taxtree_relatedness(taxtree_t* tree, int64_t id1, int64_t id2);

void taxtree_save(taxtree_t* tree, PF_t* file);
taxtree_t* taxtree_load(PF_t* file);

#endif
