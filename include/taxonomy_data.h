#ifndef DEF_TAXDATA
#define DEF_TAXDATA

#include <stdio.h>
#include <stdint.h>

typedef struct relation_table relation_table;

relation_table* rt_alloc();
relation_table* rt_copy(relation_table* target);
void rt_free(relation_table* target);
char* rt_taxrank(relation_table* target, int64_t taxid);

void rt_setmaxid(relation_table* target, int64_t maxid);
int rt_setparent(relation_table* target, int64_t id, int64_t parent);
int rt_settaxrank(relation_table* target, int64_t id, const char* name);

int rt_write(relation_table* target, FILE* f);
int rt_save(relation_table* target, const char* filename);
int rt_read(relation_table* target, FILE* f);
int rt_load(relation_table* target, const char* filename);

void rt_printlineage(relation_table* target, int64_t taxid, FILE* stream);

int64_t rt_parent(relation_table* target, int64_t taxid);
int64_t rt_parent_with_rank(relation_table* target, int64_t taxid, const char* rank);
int rt_is_parent_of(relation_table* target, int64_t node, int64_t parent_candidate);
int64_t is_good_taxrank(const char* name, int depth);
int64_t rt_goodparent(relation_table* target, int64_t taxid);
int64_t rt_lastcommonnode(relation_table* target, int64_t taxidA, int64_t taxidB);
void rt_bypass_badnodes(relation_table* target);
/*
note: global variables such as these are typically considered to be bad practice.
In this case, using them guarantees that testing for equality without strcmp is correct,
but it might  be a good idea to change how the process works in the future and remove this.
*/
extern const char* TAXREL_REPLICATES;
extern const char* TAXREL_SAME_SUBSPECIES;
extern const char* TAXREL_SAME_SPECIES;
extern const char* TAXREL_SAME_GENUS;
extern const char* TAXREL_SAME_FAMILY;
extern const char* TAXREL_SAME_ORDER;
extern const char* TAXREL_SAME_CLASS;
extern const char* TAXREL_SAME_PHYLUM;
extern const char* TAXREL_SAME_KINGDOM;
extern const char* TAXREL_DIFFERENT;

const char* rt_get_relation(relation_table* target, int64_t taxidA, int64_t taxidB, int maxsim);

typedef struct taxonomy_data taxdata;

taxdata* taxdata_alloc();
void taxdata_free(taxdata* target);

int taxdata_save(taxdata* target, const char* filename);
int taxdata_load(taxdata* target, const char* filename);

int taxdata_set_relations(taxdata* target, relation_table* rt);
int taxdata_set_name_nocheck(taxdata* target, int64_t taxid, const char* name);
char* taxdata_get_name(taxdata* target, int64_t taxid);
int64_t taxdata_get_taxid(taxdata* target, const char* name);
int taxdata_embed_relations(taxdata* target, relation_table* rt);
/*int taxdata_set_taxnames(taxdata* target, int64_t* taxids, char** taxnames, size_t amount);
int64_t taxdata_get_taxid(taxdata* target, const char* name);
relation_table* taxdata_relations(taxdata* target);*/
/**/
int taxdata_generate_subset(taxdata* target, const char* subset_code, const char* listfn, const char* taxfn, FILE* out);
int taxdata_generate_subset_interactive(taxdata* target, const char* listfn, const char* taxfn, FILE* out);
void taxdata_printlineage(taxdata* target, int64_t taxid, FILE* stream);


#endif