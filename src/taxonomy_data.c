#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "taxonomy_data.h"
#include "datamap.h"

struct relation_table {
    int64_t* parent;
    char** taxrank;
    int64_t nodecount;
};
relation_table* rt_alloc() {
    relation_table* result;
    result = malloc(sizeof(relation_table));
    result->parent = NULL;
    result->taxrank = NULL;
    result->nodecount = 0;
    return result;
}
relation_table* rt_copy(relation_table* target) {
    relation_table* result;
    int64_t i;
    if (!target)return NULL;

    result = malloc(sizeof(relation_table));
    result->parent = (int64_t*)malloc(sizeof(int64_t)*target->nodecount) ;
    result->taxrank = (char**)malloc(sizeof(char*)*target->nodecount);
    result->nodecount = target->nodecount;

    memcpy(result->parent, target->parent, sizeof(int64_t)*target->nodecount);
    for (i = 0;i < target->nodecount;i++) {
        if (target->taxrank[i]) {
            result->taxrank[i] = malloc(strlen(target->taxrank[i]) + 1);
            memcpy(result->taxrank[i], target->taxrank[i], strlen(target->taxrank[i]) + 1);
        }
        else {
            result->taxrank[i] = NULL;
        }
    }
    return result;
}
void rt_free(relation_table* target) {
    int64_t i;
    if (target->taxrank) {
        for (i = 0;i < target->nodecount;i++) {
            if (target->taxrank)
                free(target->taxrank[i]);
            target->taxrank[i] = NULL;
        }
        free(target->taxrank);
        target->taxrank = NULL;
    }
    if (target->parent) {
        free(target->parent);
        target->parent = NULL;
    }
    target->nodecount = 0;
    free(target);
}
int64_t _rt_map(relation_table* target, int64_t taxid) {
    /* This may get more complicated if the underlying data becomes too sparse for array-based storage */
    if (taxid < 0) return -taxid;
    if (taxid > target->nodecount) return 0;
    return taxid;
}
char* rt_taxrank(relation_table* target, int64_t taxid) {
    return target->taxrank[_rt_map(target, taxid)];
}

void rt_setmaxid(relation_table* target, int64_t maxid) {
    if (!target)return;
    target->nodecount = maxid+1;
    target->parent = (int64_t*)malloc(sizeof(int64_t)*target->nodecount);
    memset(target->parent, 0, sizeof(int64_t)*target->nodecount);
    target->taxrank = (char**)malloc(sizeof(char*)*target->nodecount);
    memset(target->taxrank, 0, sizeof(char*)*target->nodecount);
}
int rt_setparent(relation_table* target, int64_t id, int64_t parent) {
    int64_t curid;
    curid = _rt_map(target, id);
    if (!target || id < 0 || curid >= target->nodecount)return 1;
    target->parent[curid] = parent;
    return 0;
}
int rt_settaxrank(relation_table* target, int64_t id, const char* name) {
    size_t slen;
    int64_t curid;
    curid = _rt_map(target, id);
    if (!target || id < 0 || curid >= target->nodecount)return 1;
    slen = strlen(name);
    target->taxrank[curid] = (char*)malloc(slen + 1);
    memcpy(target->taxrank[curid], name, slen + 1);
    return 0;
}

int rt_write(relation_table* target, FILE* f) {
    int32_t i32;
    int64_t i;
    fwrite(&(target->nodecount), sizeof(int64_t), 1, f);
    for (i = 0;i < target->nodecount;i++) {
        fwrite(target->parent + i, sizeof(int64_t), 1, f);
        i32 = 0;
        if (target->taxrank[i])
            i32 = (int32_t) strlen(target->taxrank[i]);
        fwrite(&i32, sizeof(int32_t), 1, f);
        if (i32 > 0)
            fwrite(target->taxrank[i], 1, (size_t)i32, f);
    }
    return 0;
}
int rt_save(relation_table* target, const char* filename) {
    FILE* f;
    int32_t i32;
    if (!target)return 2;
#ifdef _WIN32
    fopen_s(&f, filename, "wb");
#else
    f = fopen(filename, "wb");
#endif
    if (!f)return 1;

    fwrite("TAXN", 1, 4, f);
    i32 = 1;
    fwrite(&i32, sizeof(int32_t), 1, f);
    rt_write(target, f);
    fclose(f);
    return 0;
}
int rt_read(relation_table* target, FILE* f) {
    size_t nread;
    int32_t i32;
    int64_t i;
    nread = fread(&(target->nodecount), sizeof(int64_t), 1, f);
    rt_setmaxid(target, target->nodecount-1);
    for (i = 0; i < target->nodecount; i++) {
        nread = fread(target->parent + i, sizeof(int64_t), 1, f);
        if (nread == 0) {
            return 1;
        }
        nread = fread(&i32, sizeof(int32_t), 1, f);
        if (i32 > 0) {
            target->taxrank[i] = (char*)malloc(sizeof(char)*(i32 + 1));
            nread = fread(target->taxrank[i], 1, (size_t)i32, f);
            target->taxrank[i][i32] = 0;
        }
    }
    return 0;
}
int rt_load(relation_table* target, const char* filename) {
    FILE* f;
    int32_t i32;
    size_t nread;
    int errcode;
    char* str_in;
    if (!target)return 2;
#ifdef _WIN32
    fopen_s(&f, filename, "rb");
#else
    f = fopen(filename, "rb");
#endif
    if (!f)return 1;
    errcode = 0;
    str_in = malloc(5);
    nread = fread(str_in, 1, 4, f);
    str_in[4] = 0;
    if (nread != 4 || strcmp(str_in, "TAXN"))errcode = 3;
    if (!errcode) {
        nread = fread(&i32, sizeof(int32_t), 1, f);
        if (nread != 1 || i32 != 1) errcode = 4;
    }
    if (!errcode) {
        rt_read(target, f);
    }
    fclose(f);
    return errcode;
}

void rt_printlineage(relation_table* target, int64_t taxid, FILE* stream) {
    int64_t parent;
    parent = taxid;
#ifdef _WIN32
    fprintf(stream, "%lld", taxid);
#else
    fprintf(stream, "%zd", taxid);
#endif
    parent = rt_parent(target, taxid);
    while (parent != 0) {
#ifdef _WIN32
        fprintf(stream, ",%lld", parent);
#else
        fprintf(stream, ",%zd", parent);
#endif
        parent = rt_parent(target, parent);
    }
    fputc('\n', stream);
}


int64_t rt_parent(relation_table* target, int64_t taxid) {
    if (target->nodecount == 0)return 0;
    if (taxid <= 1)return 0; /* override the loop */
    return target->parent[_rt_map(target, taxid)];
}
int64_t rt_parent_with_rank(relation_table* target, int64_t taxid, const char* rank) {
    int64_t parentcandidate;
    if (taxid > target->nodecount || taxid < 0)return 0;
    parentcandidate = rt_parent(target, taxid);
    while (parentcandidate > 1 && strcmp(rt_taxrank(target, parentcandidate), rank)!=0) {
        parentcandidate = rt_parent(target, parentcandidate);
    }
    if (parentcandidate > target->nodecount) return 0;
    return parentcandidate;
}
int rt_is_parent_of(relation_table* target, int64_t node, int64_t parent_candidate) {
    if (!target || node<1 || node >= target->nodecount)return 0;
    while ((node = rt_parent(target, node)) != 0) {
        if (node == parent_candidate) return 1;
    }
    return 0;
}
int64_t is_good_taxrank(const char* name, int depth) {
    if (!name) return 0;
    return (
        (depth >= 8 && strcmp(name, "subspecies") == 0) ||
        (depth >= 7 && strcmp(name, "species") == 0) ||
        (depth >= 6 && strcmp(name, "genus") == 0) ||
        (depth >= 5 && strcmp(name, "family") == 0) ||
        (depth >= 4 && strcmp(name, "order") == 0) ||
        (depth >= 3 && strcmp(name, "class") == 0) ||
        (depth >= 2 && strcmp(name, "phylum") == 0) ||
        (depth >= 1 && strcmp(name, "kingdom") == 0)
        );
}
int64_t _rt_goodparent_bylocation(relation_table* target, int64_t location) {
    /* note that if the target node has no proper parent, this should return 1 or 0 */
    int64_t parentcandidate;
    parentcandidate = target->parent[location];
    while (parentcandidate > 1 && !is_good_taxrank(rt_taxrank(target, parentcandidate),10)) {
        parentcandidate = rt_parent(target, parentcandidate);
    }
    return parentcandidate;
}
int64_t rt_goodparent(relation_table* target, int64_t taxid) {
    /* note that if the target node has no proper parent, this should return 1 or 0 */
    int64_t parentcandidate;
    if (taxid > target->nodecount || taxid < 0)return 0;
    parentcandidate = rt_parent(target, taxid);
    while (parentcandidate > 1 && !is_good_taxrank(rt_taxrank(target, parentcandidate),10)) {
        parentcandidate = rt_parent(target, parentcandidate);
    }
    return parentcandidate;
}
int64_t rt_lastcommonnode(relation_table* target, int64_t taxidA, int64_t taxidB) {
    int found;
    int64_t result;
    int64_t testagainst;
    result = taxidA;
    found = 0;
    while (!found && result != 0) {
        testagainst = taxidB;
        while (!found && testagainst > 0) {
            if (testagainst == result)found = 1;
            testagainst = rt_parent(target, testagainst);
        }
        if (!found) {
            result = rt_parent(target, result);
        }
    }
    return result;
}
void rt_bypass_badnodes(relation_table* target) {
    int64_t i;
    for (i = 0;i < target->nodecount;i++) {
        target->parent[i] = _rt_goodparent_bylocation(target, i);
    }
}

const char* rt_get_relation(relation_table* target, int64_t taxidA, int64_t taxidB, int maxsim) {
    int64_t lcn;
    char* taxrank;
    lcn = rt_lastcommonnode(target, taxidA, taxidB);
    while (!is_good_taxrank(rt_taxrank(target, lcn), maxsim)) {
        if (lcn <= 1) return TAXREL_DIFFERENT;
    }
    taxrank = rt_taxrank(target, lcn);
    switch (taxrank[0]) {
    case 's':
        if (taxrank[1] == 'u')return TAXREL_SAME_SUBSPECIES;
        return TAXREL_SAME_SPECIES;
    case 'g':
        return TAXREL_SAME_GENUS;
    case 'f':
        return TAXREL_SAME_FAMILY;
    case 'o':
        return TAXREL_SAME_ORDER;
    case 'c':
        return TAXREL_SAME_CLASS;
    case 'p':
        return TAXREL_SAME_PHYLUM;
    default:
        return TAXREL_DIFFERENT;
    }
}

struct taxonomy_data {
    DM64_t* id2name;
    DM64_t* name2id;
    char** names;
    int64_t nnames;
    int64_t nallocnames;
    relation_table* rt;
    DM64_t* output_name2id;
    char** output_names;
    int64_t* output_taxids;
    int64_t noutput;
};

taxdata* taxdata_alloc() {
    taxdata* result;
    result = malloc(sizeof(taxdata));

    result->name2id = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
    result->id2name = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
    result->names = NULL;
    result->nnames = 0;
    result->nallocnames = 0;

    result->rt = rt_alloc();
    
    result->output_name2id = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
    result->output_names = NULL;
    result->output_taxids = NULL;
    result->noutput = 0;

    return result;
}
void _taxdata_unlink(taxdata* target, int finalstep) {
    int64_t i;
    free_DM64(target->output_name2id);
    target->output_name2id = NULL;

    if (target->output_names) {
        for (i = 0;i < target->noutput;i++) {
            if (target->output_names[i]) free(target->output_names[i]);
            target->output_names[i] = NULL;
        }
        free(target->output_names);
        target->output_names = NULL;
    }
    if (target->output_names) {
        free(target->output_taxids);
        target->output_taxids = NULL;
    }
    if (!finalstep) {
        target->output_name2id = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
        target->noutput = 0;
    }
}
void taxdata_free(taxdata * target) {
    int64_t i;
    if (!target) return;
    free_DM64(target->name2id);
    target->name2id = NULL;
    free_DM64(target->id2name);
    target->id2name = NULL;
    if (target->names) {
        for (i = 0;i < target->nnames;i++) {
            if (target->names[i]) free(target->names[i]);
            target->names[i] = NULL;
        }
        free(target->names);
        target->names = NULL;
    }
    rt_free(target->rt);
    target->rt = NULL;
    _taxdata_unlink(target,1);
}

void taxdata_write(taxdata* target, FILE* file) {
    int64_t i;
    int32_t i32;
    DM64_write(target->id2name, file);
    DM64_write(target->name2id, file);
    fwrite(&(target->nnames), sizeof(int64_t), 1, file);
    for (i = 0;i < target->nnames;i++) {
        i32 = (int32_t) strlen(target->names[i])+1;
        fwrite(&i32, sizeof(int32_t), 1, file);
        fwrite(target->names[i], sizeof(char), (size_t)i32, file);
    }
    rt_write(target->rt, file);
    /* note:
    output_* variables are temporary in nature and not written, but in case they
    might be int he future, a placeholder value is put instead
    */
    fputc(0, file);
}
void taxdata_read(taxdata* target, FILE* file) {
    int64_t i;
    int32_t i32;
    DM64_read(target->id2name, file);
    DM64_read(target->name2id, file);

    for (i = 0;i < target->nnames;i++) {
        if (target->names[i])free(target->names[i]);
        target->names[i] = NULL;
    }
    free(target->names);
    target->names = NULL;

    fread(&(target->nnames), sizeof(int64_t), 1, file);
    target->nallocnames = target->nnames;
    target->names = (char**)calloc(target->nallocnames, sizeof(char*));
    for (i = 0;i < target->nnames;i++) {
        fread(&i32, sizeof(int32_t), 1, file);
        target->names[i] = (char*)malloc(sizeof(char)*i32);
        fread(target->names[i], sizeof(char), (size_t)i32, file);
        target->names[i][i32 - 1] = 0;
    }
    if (target->rt) rt_free(target->rt);
    target->rt = rt_alloc();
    rt_read(target->rt, file);
    /* note:
    output_* variables are temporary in nature and not written, but in case they
    might be int he future, a placeholder value is put instead
    */
    i32 = fgetc(file);
    if (i32) fprintf(stderr, "unspported");
}
int taxdata_save(taxdata* target, const char* filename) {
    FILE* f;
    int32_t i32;
    if (!target)return 2;
#ifdef _WIN32
    fopen_s(&f, filename, "wb");
#else
    f = fopen(filename, "wb");
#endif
    if (!f)return 1;

    fwrite("TAXD", 1, 4, f);
    i32 = 1;
    fwrite(&i32, sizeof(int32_t), 1, f);
    taxdata_write(target, f);
    fclose(f);
    return 0;
}
int taxdata_load(taxdata* target, const char* filename){
    FILE* f;
    int32_t i32;
    size_t nread;
    int errcode;
    char* str_in;
    if (!target)return 2;
#ifdef _WIN32
    fopen_s(&f, filename, "rb");
#else
    f = fopen(filename, "rb");
#endif
    if (!f)return 1;
    errcode = 0;
    str_in = malloc(5);
    nread = fread(str_in, 1, 4, f);
    str_in[4] = 0;
    if (nread != 4 || strcmp(str_in, "TAXD"))errcode = 3;
    if (!errcode) {
        nread = fread(&i32, sizeof(int32_t), 1, f);
        if (nread != 1 || i32 != 1) errcode = 4;
    }
    if (!errcode) {
        taxdata_read(target, f);
    }
    fclose(f);
    return errcode;
}

int taxdata_set_relations(taxdata* target, relation_table* rt) {
    if (target && rt) {
        rt_free(target->rt);
        target->rt = rt_copy(rt); /* avoid free'ing the original rt by accident */
        return 0;
    }
    return 1;
}
int taxdata_embed_relations(taxdata* target, relation_table* rt) {
    if (target && rt) {
        rt_free(target->rt);
        target->rt = rt; /* the original will be free'd when this tsruct is free'd */
        return 0;
    }
    return 1;
}
int taxdata_set_name_nocheck(taxdata* target, int64_t taxid, const char* name) {
    int64_t nid;
    size_t nlen;
    if (!target)return 1;
    nlen = strlen(name);
    nid = target->nnames;
    if (nid >= target->nallocnames) {
        if (target->nallocnames < 0x100)
            target->nallocnames = 0x100;
        else if (target->nallocnames < 0x10000000)
            target->nallocnames *= 2;
        else
            target->nallocnames += 0x10000000;
        target->names = (char**)realloc(target->names, sizeof(char*)*target->nallocnames);
    }
    target->nnames++;
    DM64_append(target->id2name, (char*)(&taxid), sizeof(int64_t), nid);
    DM64_append(target->name2id, name, (int)nlen, taxid); /* note that adding duplicate names may invalidate this datamap */
    target->names[nid] = (char*) malloc(nlen+1);
    memcpy(target->names[nid], name, nlen + 1);
    return 0;
}

char* taxdata_get_name(taxdata* target, int64_t taxid) {
    int64_t nid;
    int nullflag;
    if (!target)return NULL;
    if (taxid == 0)return NULL;
    nid = DM64_get(target->id2name, (char*)(&taxid), sizeof(int64_t), &nullflag);
    if (nullflag)
        return NULL;
    return target->names[nid];
}
int64_t taxdata_get_taxid(taxdata* target, const char* name) {
    int64_t taxid;
    int nullflag;
    if (!target)return 0;
    taxid = DM64_get(target->name2id, name , (int) strlen(name), &nullflag);
    if (nullflag)
        return 0;
    return taxid;
}
int followed_by(const char* source, const char* query) {
    int result = 0;
    size_t l;
    l = 0;
    if (!source || !query) return 0;
    while (source[l] && query[l] && query[l]==source[l]) {
        l++;
    }
    if (query[l] == 0)return 1;
    return 0;
}
int is_separating_character(char c) {
    return (c == ' ' || c == '\t' || c == ';' || c == '\n' || c=='\r');
}
size_t wordend(const char* source, size_t startat) {
    while (source[startat] != 0 && !is_separating_character(source[startat])) {
        startat++;
    }
    return startat;
}
size_t nextword(const char* source, size_t startat) {
    while (source[startat] != 0 && !is_separating_character(source[startat])) {
        startat++;
    }
    while (source[startat] != 0 && is_separating_character(source[startat])) {
        startat++;
    }
    return startat;
}

void _sel_select(uint32_t* selected, int64_t id) {
    selected[id] |= 0b001;
}
void _sel_unselect(uint32_t* selected, int64_t id) {
    if (selected[id] & 0b001)selected[id] -= 0b001;
}
void _sel_lock(uint32_t* selected, int64_t id) {
    selected[id] |= 0b010;
}
void _sel_unlock(uint32_t* selected, int64_t id) {
    if (selected[id] & 0b010)selected[id] -= 0b010;
}
void _sel_add_candidate(uint32_t* selected, int64_t id) {
    selected[id] |= 0b100;
}
void _sel_remove_candidate(uint32_t* selected, int64_t id) {
    if (selected[id] & 0b100)selected[id] -= 0b100;
}
void _sel_mark(uint32_t* selected, int64_t id) {
    selected[id] |= 0b1000;
}
void _sel_unmark(uint32_t* selected, int64_t id) {
    if (selected[id] & 0b1000)selected[id] -= 0b1000;
}

int _sel_is_selected(uint32_t* selected, int64_t id) {
    return selected[id] & 0b001;
}
int _sel_is_locked(uint32_t* selected, int64_t id) {
    return selected[id] & 0b010;
}
int _sel_is_candidate(uint32_t* selected, int64_t id) {
    return selected[id] & 0b100;
}
int _sel_is_marked(uint32_t* selected, int64_t id) {
    return selected[id] & 0b1000;
}

void _taxdata_generate_subset_random(taxdata* target, uint32_t* selected, int64_t maxcount) {
    int64_t* possibilities;
    int64_t tmpcount,i,k,tmp;
    int64_t rndres;
    tmpcount = 0;
    for (i = 0;i < target->noutput;i++) {
        if (_sel_is_candidate(selected,i)) {
            tmpcount++;
        }
    }
    possibilities = (int64_t*)malloc(sizeof(int64_t)*tmpcount);
    k = 0;
    for (i = 0;i < target->noutput;i++) {
        if (_sel_is_candidate(selected, i)) {
            possibilities[k] = i;
            k++;
        }
    }
    /*
    we want to select a random subset of n elements from the available ones
    depending on the value we will use a different approach
    */
    if (maxcount > tmpcount / 2) {
        /* randomly unselect elements*/
        k = tmpcount;
        while (k > maxcount) {
            rndres = rand() % k;
            tmp = possibilities[rndres];
            possibilities[rndres] = possibilities[k - 1];
            possibilities[k - 1] = tmp;
            k--;
        }
    }
    else {
        /* randomly select elements*/
        k = 0;
        while (k < maxcount) {
            rndres = (rand() % (int)(tmpcount - k)) + k;
            tmp = possibilities[rndres];
            possibilities[rndres] = possibilities[k];
            possibilities[k] = tmp;
            k++;
        }
    }
    for (k = maxcount;k < tmpcount;k++) {
        _sel_remove_candidate(selected, possibilities[k]);
    }
    free(possibilities);
}
void _taxdata_generate_subset_equalize(taxdata* target, uint32_t* selected, uint64_t* parent, int64_t maxcount) {
    int64_t i,j,k,tmpcount,rndres, tmp;
    int64_t* possibilities;
    uint8_t* markdone;
    int64_t matches;
    int stop;

    stop = 0;
    j = 0;
    matches = 0;
    markdone = (uint8_t*)malloc(sizeof(uint8_t)*target->noutput);
    memset(markdone, 0, sizeof(uint8_t)*target->noutput);
    while (!stop) {
        stop = 1;
        i = j;
        tmpcount = 0;
        while (i < target->noutput && stop == 1) {
            if ((selected[i] & 0b100) && !markdone[i]) {
                stop = 0;
                tmpcount = 1;
                j = i;
                matches++;
            }
            i++;
        }
        for (i = i;i < target->noutput;i++) {
            if (_sel_is_candidate(selected, i)) {
                if (parent[i] == parent[j])
                    tmpcount++;
            }
        }
        if (tmpcount < maxcount) {
            for (i = j;i < target->noutput;i++) {
                if (_sel_is_candidate(selected,i) && !markdone[i]){
                    if (parent[i] == parent[j])
                        _sel_remove_candidate(selected, i);
                }
            }
        }
        else {
            possibilities = (int64_t*)malloc(sizeof(int64_t)*tmpcount);
            k = 0;
            for (i = j;i < target->noutput;i++) {
                if (selected[i] & 0b100) {
                    if (parent[i] == parent[j]) {
                        possibilities[k] = i;
                        k++;
                    }
                }
            }
            /*
            we want to select a random subset of n elements from the available ones
            depending on the value we will use a different approach
            */
            if (maxcount > tmpcount / 2) {
                /* randomly unselect elements*/
                k = tmpcount;
                while (k > maxcount) {
                    rndres = rand() % k;
                    tmp = possibilities[rndres];
                    possibilities[rndres] = possibilities[k - 1];
                    possibilities[k - 1] = tmp;
                    k--;
                }
            }
            else {
                /* randomly select elements*/
                k = 0;
                while (k < maxcount) {
                    rndres = (rand() % (tmpcount - k)) + k;
                    tmp = possibilities[rndres];
                    possibilities[rndres] = possibilities[k];
                    possibilities[k] = tmp;
                    k++;
                }
            }
            for (k = 0;k < maxcount;k++) {
                markdone[possibilities[k]]=1;
            }
            for (k = maxcount;k < tmpcount;k++) {
                _sel_remove_candidate(selected, possibilities[k]);
                markdone[possibilities[k]] = 1;
            }
            free(possibilities);
        }
    }
    free(markdone);
}
int _sel_save_selected(taxdata* td, uint32_t* selected, const char* filename, FILE* out) {
    int64_t i;
    int donotclose;
    FILE* f;
    donotclose = 0;
    if (!td)return 2;
    if (out) {
        f = out;
        donotclose = 1;
    }
    else if (!filename || strcmp(filename, "stdout") == 0) {
        f = stdout;
        donotclose = 1;
    }
    else if (strcmp(filename, "stderr") == 0) {
        f = stderr;
        donotclose = 1;
    }
    else {
#ifdef _WIN32
        fopen_s(&f, filename, "wb");
#else
        f = fopen(filename, "wb");
#endif
    }
    if (!f) return 1;
    for (i = 0;i < td->noutput;i++) {
        if (_sel_is_selected(selected, i)) fprintf(f, "%s\n", td->output_names[i]);
    }
    if (!donotclose)
        fclose(f);
    return 0;
}
size_t _taxdata_generate_subset_parseexpr(taxdata* target, const char* subsetcode, size_t start, uint32_t* selected) {
    int64_t i,j,n,node;
    int nf;
    int64_t* tmparray;
    uint8_t* newcandidates;
    uint8_t* oldcandidates;
    size_t pos1,pos2,wordlen,end;
    char* tmpstr1;
    char* tmpstr2;
    /*
    expression format:
      <command> [arguments...];
    commands:
      @selected <expr>                  selects all values which match <expr> within the currently selected items
      @rm <expr>                        unselects all values matching <expr>
      @lock <expr>                      locks the state of values matching <expr>
      @unlock <expr>                    unlocks the state of values which match <expr>
      @mark <expr>                      marks the state of values matching <expr>
      @unmark <expr>                    unmarks the state of values matching <expr>
      @equalize <rank> <n> <expr>       ensures that exactly <n> elements are selected for each node with taxonomic rank <rank>
      @random <n> <expr>                selects a random subset of <n> values matches by <expr>, if <expr> matches less than <n> values,
                                        all values are selected
      @atleast <n> <expr>               selects all nodes which match expr if the total amount is at least <n>
      @not <expr>                       selects all values not matched by <expr>
      @all                              selects everything
      @select <expr>                    selects values which match <expr>
      @force <expr>                     select all values (including locked values) based on <expr>
      @same <rank> <expr>               selects everything which is within the same taxonpmic rank as the values matched by <expr>
                                        example: "@same genus 123456" would select all entries with the same genus as 123456
      @formarked <expr>                 repeatedly applies <expr> to select values based on marked selected values.
                                        example: "@formarked @atleast 20 @same genus;" 
      @forlocked <expr>                 repeatedly applies <expr> to select values based on lockked selected values. Trying to lock
                                        or unlock values within expr will result in undefined behaviour.
                                        example: "@forlocked @atleast 20 @same genus;"
      @count <expr>                     print the amount of values matched by <expr> to stdout
      @print <expr>                     print the values matches by <expr> to stdout
      @lineage <expression>             print the lineage of all entries matches by <expr>
      @done                             end the selected process, everything after this is ignored;
    selection format:
      +taxid                            selects all matching taxid
      {organism}                        selects all matching organism names (note: use "\}" to escape brackets)
      [name]                            selects all matching output names  (note: use "\]" to escape curly braces)
    internal formatting:
      XXX0: unselected
      XXX1: selected
      XX1X: locked
      X1XX: candidate
      1XXX: marked

    example code:
      @mark @equalize species 20 @all; @select @formarked @same genus;
    */
    if (subsetcode[start] == '@') {
        end = start;
        if (followed_by(subsetcode + start + 1, "count")) {
            end = nextword(subsetcode, start + 5);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    n++;
                }
            }
#ifdef _WIN32
            fprintf(stdout, "%lld\n", n);
#else
            fprintf(stdout, "%zd\n", n);
#endif
        }
        if (followed_by(subsetcode + start + 1, "print")) {
            end = nextword(subsetcode, start + 5);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    fprintf(stdout, "%s\n", target->output_names[i]);
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "namerank")) {
            pos1 = nextword(subsetcode, start + 8);
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;

            end = nextword(subsetcode, pos2 - 1);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i))
                    n++;
            }
            tmparray = (int64_t*)malloc(sizeof(int64_t)*n);
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    node = target->output_taxids[i];
                    tmpstr2 = rt_taxrank(target->rt, node);
                    while (node > 0 && (!tmpstr2 || strcmp(tmpstr2, tmpstr1) != 0)) {
                        node = rt_parent(target->rt, node);
                        tmpstr2 = rt_taxrank(target->rt, node);
                    }
                    if (node > 1) {
                        n = DM64_get(target->id2name, (char*)&node, sizeof(int64_t), &nf);
                        if (nf) {
                            fprintf(stdout, "%s\n", target->names[n]);
                        }
                        else {
                            fprintf(stdout, "Unknown %s\n", tmpstr1);
                        }
                    }
                    tmparray[n] = i;
                    n++;
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "all")) {
            end = nextword(subsetcode, start + 3);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                _sel_add_candidate(selected, i);
            }
        }
        if (followed_by(subsetcode + start + 1, "selected")) {
            end = nextword(subsetcode, start + 8);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected,i)) {
                    n++;
                    if (!_sel_is_selected(selected,i)) _sel_remove_candidate(selected, i);
                }
            }
            if (n == 0) {
                for (i = 0;i < target->noutput;i++) {
                    if (_sel_is_selected(selected, i)) _sel_add_candidate(selected, i);
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "add")) {
            end = nextword(subsetcode, start + 3);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    if (!(_sel_is_locked(selected, i))) _sel_select(selected, i);
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "rm")) {
            end = nextword(subsetcode, start + 2);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    if (!(_sel_is_locked(selected, i))) _sel_unselect(selected,i);
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "lock")) {
            end = nextword(subsetcode, start + 4);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    _sel_lock(selected,i);
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "unlock")) {
            end = nextword(subsetcode, start + 6);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    _sel_unlock(selected, i);
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "mark")) {
            end = nextword(subsetcode, start + 4);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    _sel_mark(selected, i);
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "unmark")) {
            end = nextword(subsetcode, start + 6);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    _sel_unmark(selected, i);
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "equalize")) {

            pos1 = nextword(subsetcode, start + 8);
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;

            pos1 = nextword(subsetcode, pos2);
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr2 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr2, subsetcode + pos1, wordlen);
            tmpstr2[wordlen] = 0;
            n = atoll(tmpstr2);
            free(tmpstr2);

            end = nextword(subsetcode, pos2);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);

            tmparray = (int64_t*)malloc(sizeof(int64_t)*target->noutput);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    tmpstr2 = rt_taxrank(target->rt, target->output_taxids[i]);
                    if (!tmpstr2 || strcmp(tmpstr2, tmpstr1) != 0)
                        tmparray[i] = rt_parent_with_rank(target->rt, target->output_taxids[i], tmpstr1);
                    else
                        tmparray[i] = target->output_taxids[i];
                }
            }
            free(tmpstr1);
            _taxdata_generate_subset_equalize(target, selected, tmparray, n);
            free(tmparray);
        }
        if (followed_by(subsetcode + start + 1, "random")) {
            pos1 = nextword(subsetcode, start + 6);
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;
            n = atoll(tmpstr1);
            free(tmpstr1);

            end = nextword(subsetcode, pos2);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            _taxdata_generate_subset_random(target, selected, n);
        }
        if (followed_by(subsetcode + start + 1, "seed")) {
            pos1 = nextword(subsetcode, start + 4);
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;
            n = atoll(tmpstr1);
            free(tmpstr1);
            srand((unsigned int)n);

            end = nextword(subsetcode, pos2);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
        }
        if (followed_by(subsetcode + start + 1, "atleast")) {
            pos1 = nextword(subsetcode, start + 7);
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;
            n = atoll(tmpstr1);
            free(tmpstr1);

            end = nextword(subsetcode, pos2);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);

            j = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    j++;
                }
            }
            if (j < n) {
                for (i = 0;i < target->noutput;i++) {
                    if (_sel_is_candidate(selected, i)) {
                        _sel_remove_candidate(selected, i);
                    }
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "forlocked")) {
            newcandidates = (uint8_t*)malloc(sizeof(uint8_t)*target->noutput);
            oldcandidates = (uint8_t*)malloc(sizeof(uint8_t)*target->noutput);
            memset(newcandidates, 0, sizeof(uint8_t)*target->noutput);
            for (i = 0;i < target->noutput;i++) {
                oldcandidates[i] = _sel_is_candidate(selected, i);
                _sel_remove_candidate(selected, i);
            }
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_locked(selected, i)) {
                    for (j = 0;j < target->noutput;j++) {
                        _sel_remove_candidate(selected, j);
                    }
                    _sel_add_candidate(selected, i);
                    end = nextword(subsetcode, start + 9);
                    if (subsetcode[end - 1] == ';') end = end - 1;
                    end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
                    for (j = 0;j < target->noutput;j++) {
                        if (_sel_is_candidate(selected, j))
                            newcandidates[j] = 1;
                    }
                }
            }
            for (i = 0;i < target->noutput;i++) {
                if (oldcandidates[i] || newcandidates[i])
                    _sel_add_candidate(selected, i);
            }
            free(newcandidates);
            free(oldcandidates);
        }
        if (followed_by(subsetcode + start + 1, "formarked")) {
            newcandidates = (uint8_t*)malloc(sizeof(uint8_t)*target->noutput);
            oldcandidates = (uint8_t*)malloc(sizeof(uint8_t)*target->noutput);
            memset(newcandidates, 0, sizeof(uint8_t)*target->noutput);
            for (i = 0;i < target->noutput;i++) {
                oldcandidates[i] = _sel_is_candidate(selected, i);
                _sel_remove_candidate(selected, i);
            }
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_marked(selected, i)) {
                    printf("%d\r", (int)(++n));
                    for (j = 0;j < target->noutput;j++) {
                        _sel_remove_candidate(selected, j);
                    }
                    _sel_add_candidate(selected, i);
                    end = nextword(subsetcode, start + 9);
                    if (subsetcode[end - 1] == ';') end = end - 1;
                    end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
                    for (j = 0;j < target->noutput;j++) {
                        if (_sel_is_candidate(selected, j))
                            newcandidates[j] = 1;
                    }
                }
            }
            for (i = 0;i < target->noutput;i++) {
                if (oldcandidates[i] || newcandidates[i])
                    _sel_add_candidate(selected, i);
            }
            free(newcandidates);
            free(oldcandidates);
        }
        if (followed_by(subsetcode + start + 1, "not")){
            end = nextword(subsetcode, start + 3);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i))
                    _sel_remove_candidate(selected, i);
                else
                    _sel_add_candidate(selected, i);
            }
        }
        if (followed_by(subsetcode + start + 1, "same")) {

            pos1 = nextword(subsetcode, start + 4);
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;

            end = nextword(subsetcode,pos2-1);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i))
                    n++;
            }
            tmparray = (int64_t*)malloc(sizeof(int64_t)*n);
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    tmparray[n] = i;
                    n++;
                }
            }
            for (i = 0;i < target->noutput;i++) {
                for (j = 0;j < n;j++) {
                    node = rt_lastcommonnode(target->rt, target->output_taxids[i], target->output_taxids[tmparray[j]]);
                    tmpstr2 = rt_taxrank(target->rt, node);
                    while (node > 0 && (!tmpstr2 || strcmp(tmpstr2, tmpstr1)!=0)) {
                        node = rt_parent(target->rt, node);
                        tmpstr2 = rt_taxrank(target->rt, node);
                    }
                    if (node > 1) {
                        _sel_add_candidate(selected, i);
                        break;
                    }
                }
            }
        }
        if (followed_by(subsetcode + start + 1, "force")) {
            end = nextword(subsetcode, start + 5);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i))
                    _sel_select(selected, i);
                else
                    _sel_unselect(selected, i);
            }
        }
        if (followed_by(subsetcode + start + 1, "done")) {
            end = nextword(subsetcode, start + 4);
            return 0;
        }
        if (followed_by(subsetcode + start + 1, "generate")) {

            pos1 = nextword(subsetcode, start + 8);
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;
            _sel_save_selected(target, selected, tmpstr1,NULL);
            end = nextword(subsetcode, pos2-1);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
        }
        if (followed_by(subsetcode + start + 1, "lineage")) {
            end = nextword(subsetcode, start + 7);
            if (subsetcode[end - 1] == ';') end = end - 1;
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
            n = 0;
            for (i = 0;i < target->noutput;i++) {
                if (_sel_is_candidate(selected, i)) {
                    taxdata_printlineage(target,target->output_taxids[i],stdout);
                }
            }
        }
        if (end <= start) {
            fprintf(stderr, "Unrecognized command: %s\n", subsetcode + start);
            while (subsetcode[start] && !(subsetcode[start] == ';' && subsetcode[start - 1] != '\\')) {
                start++;
            }
            end = start;
        }
    }
    else if (subsetcode[start] == '#') {
        /* this is a comment */
        start++;
        while (subsetcode[start] && !(subsetcode[start] == ';' && subsetcode[start-1] != '\\') && !subsetcode[start] == '\n') {
            start++;
        }
        if(subsetcode[start] == '\n')
            end = _taxdata_generate_subset_parseexpr(target, subsetcode, nextword(subsetcode,start), selected);
        else
            end = start;
    }
    else if (start>0 && subsetcode[start] && subsetcode[start] != ';'){
        if (subsetcode[start] == '+') {
            pos1 = start+1;
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;
            node = atoll(tmpstr1);
            free(tmpstr1);
        }
        else if (subsetcode[start] == '{') {
            start++;
            pos1 = start;
            while (subsetcode[start] != '\0' && subsetcode[start] != '}' && subsetcode[start - 1] != '\\')
                start++;
            pos2 = start;
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;
            node = DM64_get(target->name2id, tmpstr1, (int)strlen(tmpstr1), &nf);
            if (nf)node = 0;
            free(tmpstr1);
        }
        else if (subsetcode[start] == '[') {
            start++;
            pos1 = start;
            while (subsetcode[start] != '\0' && subsetcode[start] != ']' && subsetcode[start - 1] != '\\')
                start++;
            pos2 = start;
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;
            node = DM64_get(target->output_name2id, tmpstr1, (int) strlen(tmpstr1), &nf);
            if (nf)node = 0;
            free(tmpstr1);

        }
        else {
            start++;
            pos1 = start;
            pos2 = wordend(subsetcode, pos1);
            wordlen = pos2 - pos1;
            tmpstr1 = (char*)malloc(sizeof(char)*(wordlen + 1));
            memcpy(tmpstr1, subsetcode + pos1, wordlen);
            tmpstr1[wordlen] = 0;
            node = DM64_get(target->output_name2id, tmpstr1, (int) strlen(tmpstr1), &nf);
            if (nf) {
                node = DM64_get(target->name2id, tmpstr1, (int) strlen(tmpstr1), &nf);
                if (nf)node = 0;
            }
            free(tmpstr1);
        }
        end = nextword(subsetcode, pos2)-1;
        if (subsetcode[end-1] == ';') end = end - 1;
        end = _taxdata_generate_subset_parseexpr(target, subsetcode, end, selected);
        for (i = 0;i < target->noutput;i++) {
            if (rt_is_parent_of(target->rt, target->output_taxids[i], node) || target->output_taxids[i] == node)
                _sel_add_candidate(selected, i);
        }
    }
    else {
        end = start;
    }
    return end;
}
int _taxdata_load_available(taxdata* target, const char* listfn, const char* taxfn, char sep) {
    FILE* listofnames;
    FILE* taxofnames;
    char* line;
    size_t pos2;
    int64_t i;
    int nf;
#ifdef _WIN32
    fopen_s(&listofnames, listfn, "rb");
    fopen_s(&taxofnames, taxfn, "rb");
#else
    listofnames = fopen(listfn, "rb");
    taxofnames = fopen(taxfn, "rb");
#endif
    if (!taxofnames) {
        if (listofnames)fclose(listofnames);
        if (taxofnames)fclose(taxofnames);
        return 2;
    }

    while (line = readline(taxofnames)) {
        pos2 = 0;
        while (line[pos2] != 0 && line[pos2] != sep) {
            pos2++;
        }
        if (line[pos2] != 0) {
            line[pos2] = 0;
            pos2++;
            DM64_append(target->output_name2id, line, (int) strlen(line), atoll(line + pos2));
        }
        free(line);
    }
    if (listofnames)
        fclose(taxofnames);
    else {
        rewind(taxofnames);
        listofnames = taxofnames;
    }
    target->noutput = (int64_t) numlines(listofnames, 0);
    target->output_names = (char**)malloc(sizeof(char*)*target->noutput);
    target->output_taxids = (int64_t*)malloc(sizeof(int64_t)*target->noutput);
    memset(target->output_names, 0, sizeof(char*)*target->noutput);
    memset(target->output_taxids, 0, sizeof(int64_t)*target->noutput);
    for (i = 0;i < target->noutput;i++) {
        line = readline(listofnames);
        if (line[0] != 0) {
            target->output_taxids[i] = DM64_get(target->output_name2id, line, (int) strlen(line), &nf);
            if (nf) target->output_taxids[i] = 0;
            target->output_names[i] = line;
        }
        else
            free(line);
    }
    return 0;
}

int taxdata_generate_subset(taxdata* target, const char* subset_code, const char* listfn, const char* taxfn, FILE* out) {
    size_t start;
    size_t end;
    uint32_t* selected;
    int64_t i;
    int preloaded;
    int errcode;
    
    errcode = 0;
    if (!target)return 1;
    if (target->output_taxids == NULL) {
        errcode = _taxdata_load_available(target, listfn, taxfn,'\t');
        preloaded = 0;
    }
    else {
        preloaded = 1;
    }
    if (errcode != 0)return errcode;


    selected = (int*)calloc(target->noutput, sizeof(int));
    end = strlen(subset_code);
    start = 0;
    while (start < end) {
        start = _taxdata_generate_subset_parseexpr(target, subset_code, start, selected);
        /* 
        the above function modifies selected and always returns a values such that
        subset_code[start] is equal to '\0' or ';'. (';' is considered to be a word separator)
        */
        if (start > 0) {
            start = nextword(subset_code, start);
            for (i = 0;i < target->noutput;i++) {
                _sel_remove_candidate(selected, i);
            }
        }
        else break;
    }
    _sel_save_selected(target, selected, NULL, out);
    free(selected);
    if(!preloaded)
        _taxdata_unlink(target, 0);
    return 0;
}

int taxdata_generate_subset_interactive(taxdata* target, const char* listfn, const char* taxfn, FILE* out) {
    int stop;
    int reading;
    char *buffer;
    int64_t i;
    size_t nalloc;
    uint32_t* selected;
    int preloaded;
    int errcode;

    errcode = 0;
    if (!target)return 1;
    if (target->output_taxids == NULL) {
        errcode = _taxdata_load_available(target, listfn, taxfn, '\t');
        preloaded = 0;
    }
    else {
        preloaded = 1;
    }
    if (errcode != 0)return errcode;

    stop = 0;
    nalloc = 0x1000;
    buffer = calloc(nalloc,1);
    fprintf(stderr, "Running in interractive mode. Type '@done;' to confirm your selection and exit.\n");
    fprintf(stderr, "Expressions can span multiple lines and are only terminated by semicolons (;).\n");
    selected = (int*)calloc(target->noutput, sizeof(int));
    while (!stop) {
        reading = 1;
        i = 1;
        fprintf(stderr, "\n> ");
        while (reading) {
            buffer[i] = getchar();
            if (buffer[i] == ';' && buffer[i - 1] != '\\')reading = 0;
            i++;
            if (i >= (int64_t)nalloc) {
                nalloc *= 2;
                buffer = (char*)realloc(buffer, sizeof(char)*nalloc);
            }
        }
        buffer[i] = 0;
        i = 1;
        while (buffer[i] != 0 && is_separating_character(buffer[i]))i++;
        i = _taxdata_generate_subset_parseexpr(target, buffer, i, selected);
        /*
        the above function modifies selected and always returns a values such that
        subset_code[start] is equal to '\0' or ';'. (';' is considered to be a word separator)
        */
        if (i == 0) stop = 1;
        else {
            for (i = 0;i < target->noutput;i++) {
                _sel_remove_candidate(selected, i);
            }
        }
    }
    free(buffer);
    /* selection is not saved at the end, as it is assumed that @generate is used */
    /*_sel_save_selected(target, selected, out);*/
    free(selected);
    if (!preloaded)
        _taxdata_unlink(target, 0);
    return 0;
}
void taxdata_printlineage(taxdata* target, int64_t taxid, FILE* stream) {
    int64_t parent;
    char* name;
    parent = taxid;
    name = taxdata_get_name(target, taxid);
    if (name)
        fprintf(stream, "%s", name);
    else {
#ifdef _WIN32
        fprintf(stream, "%lld", taxid);
#else
        fprintf(stream, "%zd", taxid);
#endif
    }

    parent = rt_parent(target->rt, taxid);
    while (parent != 0) {
        name = taxdata_get_name(target, parent);
        if (name)
            fprintf(stream, ",%s", name);
        else {
#ifdef _WIN32
            fprintf(stream, ",%lld", parent);
#else
            fprintf(stream, "%,zd", parent);
#endif
        }
        parent = rt_parent(target->rt, parent);
    }
    fputc('\n', stream);
}