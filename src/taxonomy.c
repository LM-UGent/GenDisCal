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

#include "taxonomy.h"
#include "hash_table.h"
#include "textparsing.h"

#define TAXONKINGDOM    0
#define TAXONPHYLUM     1
#define TAXONCLASS      2
#define TAXONORDER      3
#define TAXONFAMILY     4
#define TAXONGENUS      5
#define TAXONSPECIES    6
#define TAXONSTRAIN     7
#define TAXONUNIQID     8

#define NMAJORTAXLVLS   9

struct lineage_t {
    int64_t majornames[NMAJORTAXLVLS]; /* kingdom, pylum, class, order, family, genus, species, strain, unique id */
    size_t nextranames;
    int64_t* minornames; /* extra taxonomic names such as 'superfamily' or 'subspecies' */
    float* minorlevels; /* levels associated  with the extra names, ranging from 0.0 to 8.0 */
};

void lineage_free(lineage_t* lineage) {
    if (lineage->nextranames > 0) {
        free(lineage->minorlevels);
        free(lineage->minornames);
    }
    free(lineage);
}
int lineage_majornamerelatedness(lineage_t* A, lineage_t* B) {
    size_t i;
    int lastcommon;
    int stop;
    stop = 0;
    lastcommon = -1;
    i = 0;
    while (!stop) {
        if (A->majornames[i] >= 0 && B->majornames[i] >= 0) {
            if (A->majornames[i] != B->majornames[i]) stop = 1;
            else lastcommon = (int)i;
        }
        i++;
        if (i > TAXONUNIQID) stop = 1;
    }
    return lastcommon; /* ignore minor names here */
}

int64_t lineage_idbytaxrank(lineage_t* lineage, double level) {
    if (level >= 0.0 && level <= 8.0) {
        if ((double)((int64_t)level) == level) {
            return lineage->majornames[(size_t)level];
        }
    }
    else if (level < 0.0) return 0;
    return -1;
}
double lineage_parentlevel(lineage_t* lineage, double level) {
    if(level >= 0.0 && level <= 8.0) {
        if ((double)((int64_t)level) == level) {
            level -= 1.0;
            while (level >= 0.0 && lineage->majornames[(size_t)level] == -1)
                level -= 1.0;
        }
        else {
            level = (double)((int64_t)level);
            while (level >= 0.0 && lineage->majornames[(size_t)level] == -1)
                level -= 1.0;
        }
        return level;
    }
    else if (level > 8.0) {
        level = 8.0;
        while (level >= 0.0 && lineage->majornames[(size_t)level] == -1)
            level -= 1.0;
        return level;
    }
    return -1.0;
}
int lineage_includes_id(lineage_t* lineage, int64_t id) {
    size_t i;
    for (i = 0;i < NMAJORTAXLVLS;i++) {
        if (id == lineage->majornames[i]) return 1;
    }
    for (i = 0;i < lineage->nextranames;i++) {
        if (id == lineage->minornames[i]) return 1;
    }
    return 0;
}

static char* taxlvlprefix[NMAJORTAXLVLS] = { "k__","p__","c__","o__","f__","g__","s__","t__","i__" };
static char* _make_taxonstring(char* str, size_t begin, size_t nchars, size_t taxlvl, size_t prevstart) {
    char* result;
    /*
    if (taxlvl == TAXONSPECIES) {
        if (begin > prevstart && memcmp(str + begin, str + prevstart, begin - prevstart - 1) == 0) {
            nchars -= begin - prevstart;
            begin += begin - prevstart;
        }
    }
    */
    result = (char*)malloc(nchars + 4);
    memcpy(result, taxlvlprefix[taxlvl], 3);
    memcpy(result + 3, str + begin, nchars);
    result[nchars + 3] = 0;
    return result;
}
static double _taxrank_str2num(char* str) {
    if ((str[0] == 'd' || str[0] == 'D') && strcmp(str + 1, "omain") == 0) return -1.0;
    if ((str[0] == 'k' || str[0] == 'K') && strcmp(str + 1, "ingdom") == 0) return 0.0;
    if ((str[0] == 'p' || str[0] == 'P') && strcmp(str + 1, "hylum") == 0) return 1.0;
    if ((str[0] == 'c' || str[0] == 'C') && strcmp(str + 1, "lass") == 0) return 2.0;
    if ((str[0] == 'o' || str[0] == 'O') && strcmp(str + 1, "rder") == 0) return 3.0;
    if ((str[0] == 'f' || str[0] == 'F') && strcmp(str + 1, "amily") == 0) return 4.0;
    if ((str[0] == 'g' || str[0] == 'G') && strcmp(str + 1, "enus") == 0) return 5.0;
    if (str[0] == 's' || str[0] == 'S') {
        if (str[1] == 'p') {
            if (str[2] == '.' || strcmp(str + 2, "ecies") == 0) return 6.0;
        }
        else if (str[1] == 'u') {
            if (strcmp(str + 2, "bspecies") == 0) return 6.5;
        }
        else if (str[1] == 's') {
            if (strcmp(str + 2, "p.") == 0) return 6.5;
        }
        else if (str[1] == 't') {
            if (str[2] == '.' || strcmp(str + 2, "rain") == 0) return 7.0;
        }
    }
    return -1.0;
}

char** PCOF_S_majornames(char* str, size_t* outlen) {
    size_t ncaps, nunder;
    size_t firstcap, firstunder, lastunder;
    size_t prevstart,i,j,k,taxlvl;
    size_t datalen;
    char** result;
    size_t lenp1, lenp2;
    i = 0;
    ncaps = 0;
    nunder = 0;
    while (str[i]) {
        if (str[i] >= 'A' && str[i] <= 'Z' && nunder<1) {
            if (ncaps == 0) firstcap = i;
            ncaps++;
        }
        if (str[i] == '_') {
            if (nunder == 0) firstunder = i;
            lastunder = i;
            nunder++;
        }
        i++;
    }
    datalen = i;
    /*
    PCOF_S cases:
    Name name   -   Genus species
    Name.name   -   Genus.species
    Name_name   -   Genus_species
    Name_Name   -   Family_Genus
    Name_Name.name   -   Genus_Genus.species
    NameName_Name   -   OrderFamily_Genus
    NameName_Name_name  -   OrderFamily_Genus_species
    NameName_Name_name_name  -   OrderFamily_Genus_species_uniqID
    NameName_Name_name_name  -   OrderFamily_Genus_species_strain_uniqID
    */
    if (nunder == 0) {
        /* Name name   -   Genus species */
        /* Name.name   -   Genus.species */
        i = 0;
        while (str[i] && str[i] != '.' && str[i] != ' ') i++;
        taxlvl = TAXONGENUS;
        j = 0;
        if (str[i]) {
            result = (char**)malloc(sizeof(char*) * 2);
            result[0] = _make_taxonstring(str, 0, i, taxlvl, 0);
            taxlvl--;
            i++;
            j++;
        }
        else {
            result = (char**)malloc(sizeof(char*) * 1);
            i = 0;
        }
        result[j] = _make_taxonstring(str, i, datalen-i, taxlvl, 0);
        *outlen = j + 1;
    }
    else if (ncaps==1){
        result = (char**)malloc(sizeof(char*) * 2);
        if (str[firstunder + 1]<'A' || str[firstunder + 1]>'Z') {
            /* Name_name   -   Genus_species */
            taxlvl = TAXONGENUS;
        }
        else if (memcmp(str,str+firstunder+1,firstunder)==0){
            /* Name_Name.name   -   Genus_Genus.species */
            /* Name_Name.name   -   Genus_Genus species */
            taxlvl = TAXONGENUS;
        }
        else {
            /* Name_Name   -   Family_Genus */
            taxlvl = TAXONFAMILY;
        }
        result[0] = _make_taxonstring(str, 0, firstunder, taxlvl, 0);
        taxlvl++;
        i = firstunder + 1;
        result[0] = _make_taxonstring(str, i, datalen-i,taxlvl, 0);
    }
    else if(ncaps<6){
        if (ncaps == 0) {
            /* If case for some reason there are no caps, consider the first letter as uppercase */
            ncaps = 1;
            firstcap = 0;
        }
        /* Family is always the first taxlvl to the left of the first underscore */
        lenp1 = ncaps;
        lenp2 = nunder;
        if (lenp2 > 4) lenp2 = 4;
        *outlen = lenp1 + lenp2;
        result = (char**)malloc(sizeof(char*)*(lenp1+lenp2));
        taxlvl = TAXONGENUS - lenp1;
        i = firstcap;
        for (j = 0;j < lenp1;j++) {
            k = i+1;
            while ((str[k]<'A' || str[k]>'Z') && str[k]!='_')k++;
            result[j] = _make_taxonstring(str, i, k - i, taxlvl, 0);
            taxlvl++;
            i = k;
        }
        i = firstunder + 1;
        k = i;
        while (str[k] != '_')k++;
        result[lenp1] = _make_taxonstring(str, i, k - i, TAXONGENUS, 0);
        if (nunder > 1) {
            prevstart = i;
            i = k + 1;
            k = i;
            while (str[k] != '_')k++;
            result[lenp1+1] = _make_taxonstring(str, i, k - i, TAXONSPECIES, prevstart);
        }
        if (nunder > 2) {
            i = k + 1;
            if (nunder == 3) k = datalen;
            else k = lastunder;
            result[lenp1+2] = _make_taxonstring(str, i, k - i, TAXONSTRAIN, 0);
        }
        if (nunder > 3) {
            i = lastunder + 1;
            result[lenp1+3] = _make_taxonstring(str, i, datalen - i, TAXONUNIQID, 0);
        }
    }
    else {
        result = NULL;
    }
    return result;
}
char** GTDB_majornames(char* str, size_t datalen, size_t* outlen) {
    char** result;
    size_t i;
    result = text_splitshort(str, ";",outlen);
    for(i = 0;i<*outlen;i++){
        if (result[i][0] == 'd') memcpy(result[i], taxlvlprefix[TAXONKINGDOM],3);
    }
    return result;
}

struct taxtree_t {
    ht64_t* name2id;
    char** id2name;
    double* id2taxlvl;
    int64_t* id2parentid;
    int64_t* id2labelid;
    int64_t nfilled;
    int64_t nalloc;
};

static void _taxtree_init(taxtree_t* tree, size_t nalloc) {
    int nf;
    tree->name2id = ht64_alloc_size(nalloc);
    tree->id2name = (char**)malloc(nalloc * sizeof(char*));
    tree->id2taxlvl = (double*)malloc(nalloc * sizeof(double));
    tree->id2parentid = (int64_t*)malloc(nalloc * sizeof(int64_t));
    tree->id2labelid = (int64_t*)malloc(nalloc * sizeof(int64_t));
    tree->id2name[0] = (char*)malloc(5);
    tree->id2parentid[0] = 0;
    tree->id2taxlvl[0] = -2.0;
    tree->id2labelid[0] = 0;
    memcpy(tree->id2name[0], "root", 5);
    ht64_set(tree->name2id, "root", 4, 0, &nf);
    tree->nfilled = 1;
    tree->nalloc = nalloc;
}

static inline size_t _taxtree_reserve_array_index(taxtree_t* tree, size_t nalloc, size_t newmaxindex) {
    if (newmaxindex >= nalloc) {
        if (nalloc * 2 <= newmaxindex) {
            nalloc = newmaxindex + 1;
        }
        if (nalloc < 0x100000)
            nalloc *= 2;
        else
            nalloc += 0x100000;
        tree->id2name = (char**)realloc(tree->id2name, nalloc * sizeof(char*));
        tree->id2taxlvl = (double*)realloc(tree->id2taxlvl, nalloc * sizeof(double));
        tree->id2parentid = (int64_t*)realloc(tree->id2parentid, nalloc * sizeof(int64_t));
        tree->id2labelid = (int64_t*)realloc(tree->id2labelid, nalloc * sizeof(int64_t));
    }
    return nalloc;
}

taxtree_t* taxtree_alloc() {
    taxtree_t* tree;
    size_t nalloc;
    tree = (taxtree_t*)malloc(sizeof(taxtree_t));
    nalloc = 4096;
    _taxtree_init(tree, nalloc);
    return tree;
}
taxtree_t* taxtree_fromlineages(PF_t* file, int fmt) {
    taxtree_t* tree;
    size_t taxstart;
    size_t nfilled, nalloc;
    int64_t taxid;
    int64_t previd;
    size_t i, j, n;
    char sep;
    char* line;
    char* tmpstr;
    char* name;
    char** taxonames;
    char* lineagestring;
    int nf, prevnf;
    tree = (taxtree_t*) malloc(sizeof(taxtree_t));
    nalloc = 4096;
    _taxtree_init(tree, nalloc);
    nfilled = tree->nfilled;
    if (fmt == TAXFMT_PCOF_S)sep = ',';
    else sep = '\t';
    while (line = PFreadline(file)) {
        taxstart = text_findchar(line, sep);
        name = line;
        if (taxstart >= strlen(line)) continue;
        lineagestring = line + taxstart + 1;
        line[taxstart] = 0;
        if (fmt == TAXFMT_PCOF_S)
            taxonames = PCOF_S_majornames(lineagestring, &n);
        else if (fmt == TAXFMT_GTDB)
            taxonames = GTDB_majornames(lineagestring, strlen(lineagestring), &n);
        else {
            taxonames = NULL;
            n = 0;
        }
        previd = 0;
        nf = 0;
        prevnf = 0;
        for (i = 0;i < n;i++) {
            taxid = ht64_get(tree->name2id, taxonames[i], strlen(taxonames[i]), &nf);
            /*
            A lineage should always go from highest level (e.g. domain) to lowest level (e.g. species),
            so one we reach a novel taxon, all following names should be novel as well. Adding them
            will mess up some things, but the taxon ids will remain consistent
            */
            prevnf |= nf;
            if (nf || prevnf) {
                nalloc = _taxtree_reserve_array_index(tree, nalloc, nfilled+1);
                taxid = nfilled;
                ht64_set(tree->name2id, taxonames[i], strlen(taxonames[i]), (int64_t)nfilled, &nf);
                tree->id2name[nfilled] = taxonames[i];
                tree->id2labelid[nfilled] = nfilled;
                for (j = 0; j < NMAJORTAXLVLS; j++) {
                    if (memcmp(taxonames[i], taxlvlprefix[j], 3) == 0) {
                        tree->id2taxlvl[nfilled] = (double)j;
                        break;
                    }
                }
                tree->id2parentid[nfilled] = previd;
                nfilled++;
                if (memcmp(taxonames[i] + 1, "__", 2) == 0) {
                    j = strlen(taxonames[i]);
                    tmpstr = malloc(j + 1);
                    memcpy(tmpstr, taxonames[i], j + 1);
                    tmpstr[0] = '#';
                    tmpstr[1] = '&';
                    tmpstr[2] = '@';
                    ht64_set(tree->name2id, tmpstr, j, (int64_t)nfilled, &nf);
                    tree->id2name[nfilled] = tmpstr;
                    tree->id2labelid[nfilled-1] = nfilled;
                    tree->id2parentid[nfilled] = nfilled-1;
                    nfilled++;
                }
            }
            previd = taxid;
        }
        taxid = ht64_get(tree->name2id, name, taxstart , &nf);
        if (nf) {
            nalloc = _taxtree_reserve_array_index(tree, nalloc, nfilled);
            ht64_set(tree->name2id, name, strlen(name), (int64_t)nfilled, &nf);
            tree->id2name[nfilled] = (char*)malloc(taxstart+1);
            memcpy(tree->id2name[nfilled], name, taxstart + 1);
            tree->id2parentid[nfilled] = previd;
            tree->id2labelid[nfilled] = nfilled;
            tree->id2taxlvl[nfilled] = TAXONUNIQID;
            nfilled++;
        }
        previd = taxid;
        free(line);
    }
    tree->nfilled = nfilled;
    tree->nalloc = nalloc;
    return tree;
}
taxtree_t* taxtree_fromrelations(PF_t* file, int fmt) {
    taxtree_t* tree;
    size_t nfilled, nalloc;
    size_t i;
    char* line;
    char** parts;
    char* ownid;
    char* parentid;
    char* taxrank;
    int nf;
    int64_t realparentid;
    int64_t realownid;
    size_t outcount;
    tree = (taxtree_t*)malloc(sizeof(taxtree_t));
    nalloc = 4096;
    _taxtree_init(tree, nalloc);
    nfilled = tree->nfilled;
    while (line = PFreadline(file)) {
        /* each line adds at most two ids (node+parent) */
        nalloc = _taxtree_reserve_array_index(tree, nalloc, nfilled+1);
        if (fmt == TAXFMT_NCBI) {
            parts = text_splitshort(line, "\t|\t", &outcount);
            ownid = parts[0];
            parentid = parts[1];
            taxrank = parts[2];
            for (i = 3;i < outcount;i++) {
                free(parts[i]);
            }
            realparentid = ht64_get(tree->name2id, parentid, strlen(parentid), &nf);
            if (nf) {
                ht64_set(tree->name2id, parentid, strlen(parentid), nfilled, &nf);
                realparentid = nfilled;
                tree->id2name[realparentid] = parentid;
                tree->id2labelid[realparentid] = realparentid;
                nfilled++;
            }
            else free(parentid);
            realownid = ht64_get(tree->name2id, ownid, strlen(ownid), &nf);
            if (nf) {
                ht64_set(tree->name2id, ownid, strlen(ownid), nfilled, &nf);
                realownid = nfilled;
                tree->id2name[realownid] = ownid;
                nfilled++;
            }
            else {
                tree->id2name[realownid] = NULL;
                free(ownid);
            }
            tree->id2parentid[realownid] = realparentid;
            tree->id2labelid[realownid] = realownid;
            tree->id2taxlvl[realownid] = _taxrank_str2num(taxrank);
        }
        free(line);
    }
    tree->nfilled = nfilled;
    tree->nalloc = nalloc;
    return tree;
}
void taxtree_addleaves(taxtree_t* tree, PF_t* file, int fmt) {
    char* line;
    char* name;
    char* parentid;
    char** parts;
    size_t sep;
    size_t nfilled, nalloc, tmpsz, nparts;
    size_t i;
    int64_t realparentid;
    int nf;
    nfilled = (size_t)tree->nfilled;
    nalloc = (size_t)tree->nalloc;
    while (line = PFreadline(file)) {
        nalloc = _taxtree_reserve_array_index(tree, nalloc, nfilled);
        if (fmt == TAXFMT_CSV) sep = text_findchar(line, ',');
        else if (fmt == TAXFMT_TSV)sep = text_findchar(line, '\t');
        
        if (fmt == TAXFMT_CSV || fmt == TAXFMT_TSV) {
            /* basic formats for adding samples */
            line[sep] = 0;
            parentid = line + sep + 1;
            name = line;
            parts = NULL;
            tmpsz = strlen(name);
            ht64_set(tree->name2id, name, tmpsz, nfilled, &nf);
            tree->id2name[nfilled] = malloc(tmpsz + 1);
            tree->id2taxlvl[nfilled] = TAXONUNIQID;
            tree->id2labelid[nfilled] = nfilled;
            memcpy(tree->id2name[nfilled], name, tmpsz + 1);
            realparentid = ht64_get(tree->name2id, parentid, strlen(parentid), &nf);
            if (!nf) tree->id2parentid[nfilled] = realparentid;
            else tree->id2parentid[nfilled] = 0;
            nfilled++;
        }
        else if (fmt == TAXFMT_NCBI) {
            /* adding labels */
            parts = text_splitshort(line, "\t|\t", &nparts);
            if (nparts>=4 && strcmp("scientific name", parts[3]) == 0) {
                tmpsz = strlen(parts[1]) + 3;
                name = malloc(tmpsz + 1);
                memcpy(name + 3, parts[1], tmpsz-3);
                name[0] = '#';
                name[1] = '&';
                name[2] = '@';
                /*note: #&@ (read: "identifier and at") should be interpreted as "an alternative id is" */
                name[tmpsz] = 0;
                parentid = parts[0];
                realparentid = ht64_get(tree->name2id, parentid, strlen(parentid), &nf);
                if (!nf) {
                    tree->id2name[nfilled] = name;
                    tree->id2parentid[nfilled] = realparentid;
                    tree->id2labelid[realparentid] = nfilled;
                    tree->id2taxlvl[nfilled] = tree->id2taxlvl[realparentid]+0.001;
                    ht64_set(tree->name2id, name, tmpsz, nfilled, &nf);
                    nfilled++;
                }
                else {
                    free(name);
                }
            }
            if (parts) {
                for (i = 0;i < nparts;i++)free(parts[i]);
                free(parts);
            }
        }
        free(line);
    }
    tree->nfilled = nfilled;
    tree->nalloc = nalloc;
}

void taxtree_rename(taxtree_t* tree, PF_t* file, int fmt) {
    char* line;
    char* parts[2];
    size_t linelen;
    size_t namelen;
    size_t newnamelen;
    size_t i;
    int64_t id;
    int64_t newid;
    int nf1, nf2;
    while (line = PFreadline(file)) {
        /* old name - new name */
        parts[0] = line;
        linelen = strlen(line);
        switch (fmt) {
        case TAXFMT_CSV:
            i = text_findchar(line, ',');
            namelen = i;
            i++;
            while (line[i] != ',' && line[i]!=0 )i++;
            newnamelen = i - namelen - 1;
            break;
        case TAXFMT_TSV:
            i = text_findchar(line, '\t');
            namelen = i;
            i++;
            while (line[i] != '\t' && line[i] != 0)i++;
            newnamelen = i - namelen - 1;
            break;
        case TAXFMT_NCBI:
            i = 2;
            while (i < linelen && (line[i - 2] != '\t' || line[i - 1] != '|' || line[i] != '\t')) {
                i += text_findchar(line + i+1, '\t')+1;
            }
            namelen = i-2;
            i++;
            while (line[i] != '\t' && line[i] != '\r' && line[i] != '\n' && line[i] != 0)i++;
            newnamelen = i - namelen - 3;
            break;
        default:
            i = linelen;
            namelen = linelen;
            break;
        }
        line[namelen] = 0;
        i++;
        if (namelen < linelen) {
            parts[1] = line + namelen + 1;
            parts[1][newnamelen] = 0;
            id = ht64_get(tree->name2id, parts[0], namelen, &nf1);
            newid = ht64_get(tree->name2id, parts[1], newnamelen, &nf2);
            /* The old name can still be used - but will redirect to the new entry */
            if (nf1==0 && nf2==0) {
                ht64_set(tree->name2id, parts[0], namelen, newid, &nf1);
                free(tree->id2name[id]);
                tree->id2name[id] = malloc(newnamelen + 1);
                memcpy(tree->id2name[id], tree->id2name[newid], newnamelen + 1);
                tree->id2taxlvl[id] = tree->id2taxlvl[newid];
                tree->id2parentid[id] = tree->id2parentid[newid];
                tree->id2labelid[id] = tree->id2labelid[newid];
            }
            else if (nf2 && !nf1) {
                /* the 'new' name was not in the database - just change the name of the old one */
                ht64_set(tree->name2id, parts[0], namelen, id, &nf1);
                ht64_set(tree->name2id, parts[1], newnamelen, id, &nf1);
                free(tree->id2name[id]);
                tree->id2name[id] = malloc(newnamelen + 1);
                memcpy(tree->id2name[id], parts[1], newnamelen + 1);
            }
            else if (nf1 && !nf2) {
                /* the old name wasn't in the database before - just redirect it to the new one */
                ht64_set(tree->name2id, parts[0], namelen, newid, &nf1);
            }
            else {
                /* neither the old name nor the new name were in the database bfore - create a dummy entry and let both names redirect to it*/
                tree->nalloc = _taxtree_reserve_array_index(tree, tree->nalloc, tree->nfilled);
                newid = tree->nfilled;
                tree->id2labelid[newid] = newid;
                tree->id2name[newid] = malloc(newnamelen + 1);
                memcpy(tree->id2name[newid], parts[1], newnamelen + 1);
                tree->id2taxlvl[newid] = -1.0;
                tree->id2parentid[newid] = 0;
                tree->nfilled++;
                ht64_set(tree->name2id, parts[0], namelen, newid, &nf1);
                ht64_set(tree->name2id, parts[1], newnamelen, newid, &nf2);
            }
        }
        free(line);
    }
}
void taxtree_free(taxtree_t* target) {
    int64_t i;
    for (i = 0;i < target->nfilled;i++) {
        if (target->id2name[i])free(target->id2name[i]);
        target->id2name[i] = NULL;
    }
    free(target->id2parentid);
    free(target->id2taxlvl);
    free(target->id2name);
    ht64_free(target->name2id);
    free(target);
}
int64_t taxtree_getid(taxtree_t* tree, char* name) {
    int nf;
    int64_t result;
    size_t tmplen;
    char* tmp;
    result = ht64_get(tree->name2id,name,strlen(name),&nf);
    if (nf) {
        tmplen = strlen(name);
        tmp = malloc(tmplen + 4);
        tmp[0] = '#';
        tmp[1] = '&';
        tmp[2] = '@';
        memcpy(tmp + 3, name, tmplen);
        tmp[tmplen + 3] = 0;
        result = ht64_get(tree->name2id, tmp, tmplen+3, &nf);
        free(tmp);
        if (nf)
            return -1;
        else
            result = tree->id2parentid[result];
    }
    return result;
}
static void _taxtree_fill_lineage_by_id(lineage_t* target, taxtree_t* tree, int64_t id) {
    int64_t curid;
    size_t i;
    double linlvl;
    curid = id;
    for (i = 1;i < NMAJORTAXLVLS;i++) {
        target->majornames[i] = -1;
    }
    target->nextranames = 0;
    while (curid > 0) {
        linlvl = tree->id2taxlvl[curid];
        /* check if it is an integer*/
        if ((double)((size_t)linlvl) == linlvl && linlvl>=0) {
            target->majornames[(size_t)linlvl] = curid;
        }
        if (curid == tree->id2parentid[curid]) curid = 0;
        else curid = tree->id2parentid[curid];
    }
}
int64_t taxtree_getparentid(taxtree_t* tree, int64_t id) {
    return tree->id2parentid[id];
}
lineage_t* taxtree_lineageof(taxtree_t* tree, int64_t id) {
    lineage_t* result;
    result = NULL;
    if (id >= 0 && id < tree->nfilled) {
        result = calloc(1, sizeof(lineage_t));
        _taxtree_fill_lineage_by_id(result, tree, id);
    }
    return result;
}
char* taxtree_getnameptr(taxtree_t* tree, int64_t id) {
    if (id >= 0 && id < tree->nfilled)
        return tree->id2name[id];
    else
        return "?";
}
char* taxtree_getlabelptr(taxtree_t* tree, int64_t id) {
    int nf = 0;
    if (id >= 0 && id < tree->nfilled) {
        if (id > 0 && tree->id2name[tree->id2labelid[id]] && memcmp(tree->id2name[tree->id2labelid[id]], "#&@", 3) == 0)
            return tree->id2name[tree->id2labelid[id]] + 3;
        else if (id > 0 && tree->id2name[id])
            return tree->id2name[ht64_get(tree->name2id, tree->id2name[id], strlen(tree->id2name[id]), &nf)];
        else
            return "organism";
    }
    else
        return NULL;
}
char* taxtree_generatereadablename(taxtree_t* tree, int64_t id) {
    int nf = 0;
    char* result;
    char* tmp[3] = { "unclassified",NULL,NULL };
    int64_t tmpid;
    size_t tokeep;
    size_t namelen;
    result = NULL;
    if (id >= 0 && id < tree->nfilled) {
        if (tree->id2taxlvl[id] < TAXONSPECIES) {
            tmp[0] = taxtree_getlabelptr(tree, id);
            tmpid = 0;
            tokeep = 1;
            if (tree->id2taxlvl[id] > TAXONPHYLUM || tree->id2taxlvl[tmpid] <= 0.0) {
                tmpid = tree->id2parentid[id];
                while (tmpid > 0 && (tree->id2taxlvl[tmpid] != TAXONPHYLUM || tree->id2taxlvl[tmpid] <= 0.0)) {
                    tmpid = tree->id2parentid[tmpid];
                }
            }
            if (tree->id2taxlvl[tmpid] == TAXONPHYLUM) {
                tmp[1] = taxtree_getlabelptr(tree,tmpid);
                tokeep = 2;
            }
            if (tokeep == 2)
                result = text_join(tmp, " (", NULL, ")", tokeep);
            else
                result = text_join(tmp, " ", NULL, NULL, tokeep);
        }
        else {
            tokeep = 0;
            tmpid = id;
            while (tmpid > 0 && (tree->id2taxlvl[tmpid] > TAXONSPECIES || tree->id2taxlvl[tmpid]<=0.0)) {
                tmpid = tree->id2parentid[tmpid];
            }
            if (tree->id2taxlvl[tmpid] == TAXONSPECIES) {
                tmp[1] = taxtree_getlabelptr(tree, tmpid);
                while (tmpid > 0 && (tree->id2taxlvl[tmpid] > TAXONGENUS || tree->id2taxlvl[tmpid] <= 0.0)) {
                    tmpid = tree->id2parentid[tmpid];
                    if (tmpid == tree->id2parentid[tmpid]) break;
                }
                tokeep = 1;
                if (tree->id2taxlvl[tmpid] == TAXONGENUS) tmp[0] = taxtree_getlabelptr(tree, tmpid);
                else {
                    tmp[0] = tmp[1];
                    tokeep = 0;
                }
                tokeep++;
                tmpid = id;
                while (tree->id2parentid[tmpid] > 0 && tree->id2taxlvl[tree->id2parentid[tmpid]] > TAXONSPECIES) {
                    tmpid = tree->id2parentid[tmpid];
                }
                tmp[tokeep] = taxtree_getlabelptr(tree, tmpid);
                tokeep++;
            }
            else {
                tokeep = 2;
                tmp[1] = taxtree_getlabelptr(tree, tmpid);
            }
            for (tmpid = tokeep - 1;tmpid > 0;tmpid--) {
                namelen = strlen(tmp[tmpid - 1]);
                if (memcmp(tmp[tmpid], tmp[tmpid-1], namelen) == 0) {
                    if (iswordsep(tmp[tmpid][namelen]))tmp[tmpid] += namelen+1;
                }
                else if (memcmp(tmp[tmpid] + 1, "__", 2) == 0) {
                    if (memcmp(tmp[tmpid]+1, tmp[tmpid - 1]+1, namelen-1) == 0) {
                        if (iswordsep(tmp[tmpid][namelen]))tmp[tmpid] += namelen + 1;
                    }
                }
            }
            for (tmpid = 0;tmpid < (int64_t)tokeep;tmpid++) {
                if (memcmp(tmp[tmpid] + 1, "__", 2) == 0) {
                    if (tmp[tmpid])tmp[tmpid] += 3;
                }
            }
            result = text_join(tmp, " ", NULL, NULL, tokeep);
        }
    }
    return result;
}
int taxtree_relatedness(taxtree_t* tree, int64_t id1, int64_t id2)
{
    lineage_t lin1;
    lineage_t lin2;
    _taxtree_fill_lineage_by_id(&lin1, tree, id1);
    _taxtree_fill_lineage_by_id(&lin2, tree, id2);
    return lineage_majornamerelatedness(&lin1,&lin2);
}

void taxtree_save(taxtree_t* tree, PF_t* file) {
    int64_t i;
    int64_t namelen;
    PFputint64(file, tree->nfilled);
    ht64_save(tree->name2id, file);
    for (i = 0;i < tree->nfilled;i++) {
        PFputint64(file, tree->id2parentid[i]);
        PFputint64(file, tree->id2labelid[i]);
        PFputdouble(file, tree->id2taxlvl[i]);
        if (tree->id2name[i]) {
            namelen = strlen(tree->id2name[i]);
            PFputint64(file, namelen);
            PFwrite(tree->id2name[i], 1, namelen, file);
        }
        else {
            PFputint64(file, 0);
        }
    }
}
taxtree_t* taxtree_load(PF_t* file) {
    taxtree_t* tree;
    size_t i;
    size_t namelen;
    tree = (taxtree_t*)malloc(sizeof(taxtree_t));
    tree->nfilled = PFgetint64(file);
    tree->nalloc = tree->nfilled;
    tree->name2id = ht64_load(file);
    tree->id2name = (char**) malloc(sizeof(char*)*tree->nalloc);
    tree->id2parentid = (int64_t*)malloc(sizeof(int64_t)*tree->nalloc);
    tree->id2taxlvl = (double*)malloc(sizeof(double)*tree->nalloc);
    tree->id2labelid = (int64_t*)malloc(sizeof(int64_t)*tree->nalloc);
    for (i = 0;i < (size_t)(tree->nfilled);i++) {
        tree->id2parentid[i] = PFgetint64(file);
        tree->id2labelid[i] = PFgetint64(file);
        tree->id2taxlvl[i] = PFgetdouble(file);
        namelen = PFgetint64(file);
        if (namelen > 0) {
            if (namelen > 10000) {
                tree->nfilled = i-1;
                fprintf(stderr, "Exceedingly long key detected. Tree loading will be aborted\n");
                return tree;
            }
            tree->id2name[i] = (char*)malloc(namelen+1);
            PFread(tree->id2name[i], 1, namelen, file);
            tree->id2name[i][namelen] = 0;
        }
        else {
            tree->id2name[i] = NULL;
        }
    }
    return tree;
}