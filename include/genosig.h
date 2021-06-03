#ifndef DEF_GENOSIG
#define DEF_GENOSIG

#include "taxonomy.h"
#include "osportstd.h"
#include "nucseq.h"

typedef struct genosig_t genosig_t;

genosig_t* genosig_alloc();
void genosig_cleargenomedata(genosig_t* target);
void genosig_clearsig(genosig_t* target);
void genosig_free(genosig_t* target);

size_t genosig_info_sigsize(genosig_t* target);
size_t genosig_info_length(genosig_t* target);
char* genosig_info_name(genosig_t* target);
int genosig_info_ismuxed(genosig_t* target);
int genosig_matchesname(genosig_t* target, int64_t taxid, char* literalname);
/* import/export data */
#define COPYLVL_NOTHING     0   /* nothing is copied - the passed array should be free'd outside of the structure-specific free */
#define COPYLVL_INTEGRATE   1   /* the passed array is integrated into the structure - it should not be free'd outside of the structure-specific free */
#define COPYLVL_FULL        2   /* the passed array is fully copied into the structure - it should still be free'd outside of the structure specific free */
genosig_t* genosig_fullgenome(nucseq** data, size_t nseqs, int is_single_strand, int copylvl);
genosig_t* genosig_suffixtrees(genosig_t* sig);
genosig_t* genosig_linkfilename(genosig_t* sig, char* fn);
genosig_t* genosig_genomefileref(char* path);
genosig_t* genosig_linkname(genosig_t* sig, char* name);
genosig_t* genosig_autoname(genosig_t* sig, char* filepath);
#define RENAMEBY_NODENAME   0
#define RENAMEBY_NODELABEL  1
#define RENAMEBY_SPECIES    2
genosig_t* genosig_renameby(genosig_t* sig, taxtree_t* tree, int copylvl, int mode);
genosig_t* genosig_integratelineage(genosig_t* sig, lineage_t* lineage);
genosig_t* genosig_importlineage(genosig_t* sig, taxtree_t* tree, int* nullflag);
genosig_t* genosig_makewindows(genosig_t* target, size_t winsize, size_t winskip);
genosig_t* genosig_findgenes(genosig_t* data, int is_single_strand, genosig_t* genesigs, int detectionmethod);
genosig_t* genosig_loadbin(PF_t* file);
genosig_t* genosig_import(char* filename);
char* genosig_astxt(genosig_t* sig, size_t* outsize);
char* genosig_asbin(genosig_t* sig, size_t* outsize);
void genosig_savetxt(genosig_t* sig, PF_t* file);
void genosig_savebin(genosig_t* sig, PF_t* file);
void genosig_export(genosig_t* sig, const char* filename);
genosig_t* genosig_multiplex(genosig_t** sigs, size_t count);
genosig_t** genosig_demux(genosig_t* muxedsigs, size_t* count);

/* basis */
typedef genosig_t* (*genosig_bf)(genosig_t*, size_t);
genosig_t* genosig_keepfilenameonly(genosig_t* sig, size_t unused);
genosig_t* genosig_GC(genosig_t* sig, size_t unused);
genosig_t* genosig_length(genosig_t* sig, size_t unused);
genosig_t* genosig_kmercount(genosig_t* sig, size_t k);
genosig_t* genosig_minhash(genosig_t* sig, size_t numhashes);
genosig_t* genosig_kmerfreq(genosig_t* sig, size_t k);
genosig_t* genosig_mmz_inf(genosig_t* sig, size_t k);
genosig_t* genosig_mmz(genosig_t* sig, size_t k);
genosig_t* genosig_karlinsig(genosig_t* sig, size_t k);
genosig_t* genosig_karlinsigL(genosig_t* sig, size_t k);
genosig_t* genosig_presenceabsence(genosig_t* sig, size_t unused);
genosig_t* genosig_bactspeciessig(genosig_t* sig, size_t unused);
genosig_t* genosig_bacttypesig(genosig_t* sig, size_t unused);
genosig_t* genosig_avgsig(genosig_t* sig, size_t is_weighted);
genosig_t* genosig_varsig(genosig_t* sig, size_t is_weighted);
genosig_t* genosig_keepnonredundant(genosig_t* sig, size_t k);

/* distance */
typedef double (*genosig_df)(genosig_t*, genosig_t*, any_t);
#define GENODIST_NONORM             0
#define GENODIST_NORMELTSBYDIMS     1
#define GENODIST_NORMELTSBYINVDIMS  2
#define GENODIST_NORMELTSBYAVG      3
#define GENODIST_NORMBYSEQ          4
#define GENODIST_NORMBYLEN          5
double genodist_bytaxonomy(genosig_t* A, genosig_t* B, any_t flags);
double genodist_manhattan(genosig_t* A, genosig_t* B, any_t flags);
double genodist_euclidian(genosig_t* A, genosig_t* B, any_t flags);
double genodist_pearscorr(genosig_t* A, genosig_t* B, any_t flags);
double genodist_pearscorr_unbound(genosig_t* A, genosig_t* B, any_t flags);
double genodist_rankcorr(genosig_t* A, genosig_t* B, any_t flags);
double genodist_rankcorr_unbound(genosig_t* A, genosig_t* B, any_t flags);
double genodist_jaccard(genosig_t* A, genosig_t* B, any_t flags);
double genodist_hamming(genosig_t* A, genosig_t* B, any_t flags);
double genodist_approxANI(genosig_t* A, genosig_t* B, any_t flags);
double genodist_approxANI_unbound(genosig_t* A, genosig_t* B, any_t flags);
double genodist_SVC(genosig_t* A, genosig_t* B, any_t maxdelta);
double genodist_satman(genosig_t* A, genosig_t* B, any_t maxdelta);
double genodist_PaSiTL(genosig_t* A, genosig_t* B, any_t maxdelta);
double genodist_sateucl(genosig_t* A, genosig_t* B, any_t maxdelta);
double genodist_ANI(genosig_t* A, genosig_t* B, any_t flags);
double genodist_samespecies(genosig_t* A, genosig_t* B, any_t flags);
double genodist_sametype(genosig_t* A, genosig_t* B, any_t flags);

double genodist_externANIb(genosig_t* A, genosig_t* B, any_t blastpath);
double genodist_externANIb_unbound(genosig_t* A, genosig_t* B, any_t blastpath);
double genodist_externdist(genosig_t* A, genosig_t* B, any_t progstring);

#endif

