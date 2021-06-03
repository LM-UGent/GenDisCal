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

#include "textparsing.h"
#include "hash_table.h"
#include "genosig.h"
#include "vecops.h"
#include "suffixtree.h"

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#define _SIGTYPE_NOSIG  0
#define _SIGTYPE_BYTE   1
#define _SIGTYPE_I32    2
#define _SIGTYPE_I64    3
#define _SIGTYPE_DBL    4
#define _SIGTYPE_SZT    5
#define _SIGTYPE_SUFTREEARRAY    6
#define _SIGTYPE_HASHTBL    7
#define _SIGTYPE_MUXED  8
#define _SIGTYPE_MAX    9

#define _SIGNAME_UNKNOWN    -1
#define _SIGNAME_NONE       0
#define _SIGNAME_SEQUENCE   1
#define _SIGNAME_COUNTS     2
#define _SIGNAME_FREQPROF   3
#define _SIGNAME_MMZ        4
#define _SIGNAME_KARLIN     5
#define _SIGNAME_PRESABS    6
#define _SIGNAME_HASHLIST   7
#define _SIGNAME_SUFTREE    8
#define _SIGNAME_MUXED      9
#define _SIGNAME_GC         10
#define _SIGNAME_GENOMELEN  11
#define _SIGNAME_FILENAME   12
#define _SIGNAME_KARLINL    13
#define _SIGNAME_MULTINORM  14
#define _SIGNAME_MAX        15

struct genosig_t {
    nucseq** contigs;
    size_t ncontigs;
    void* sigdata;
    size_t sigsize;
    size_t siglen;
    int16_t sigtype;
    int16_t signame;
    int sigmux;
    size_t uid;
    char* name;
    char* filepath;
    lineage_t* lineage;
    size_t window;
    size_t windowskip;
    double weight;
    int linkedgenome;
    int lock;
    int allocatedname;
    int singlestrand;
};

void genosig_cleargenomedata(genosig_t* target) {
    size_t i;
    if (!(target->linkedgenome)) {
        for (i = 0;i < target->ncontigs;i++) {
            clear_nucseq(target->contigs[i]);
            free(target->contigs[i]);
        }
        free(target->contigs);
        target->linkedgenome = 0;
    }
    target->contigs = NULL;
    target->ncontigs = 0;
}
void genosig_clearsig(genosig_t* sig) {
    size_t i;
    if (sig->sigdata) {
        if (sig->sigtype == _SIGTYPE_HASHTBL) ht64_free((ht64_t*)(sig->sigdata));
        else if (sig->sigtype == _SIGTYPE_SUFTREEARRAY) {
            for (i = 0;i < sig->siglen;i++) {
                if (((suftree_t**)(sig->sigdata))[i])
                    suftree_free(((suftree_t**)(sig->sigdata))[i]);
            }
            free(sig->sigdata);
        }
        else if (sig->sigtype == _SIGTYPE_MUXED) {
            for (i = 0;i < sig->siglen;i++) {
                genosig_free(((genosig_t**)(sig->sigdata))[i]);
            }
            free(sig->sigdata);
        }
        else free(sig->sigdata);
        sig->sigdata = NULL;
    }
}

genosig_t* genosig_alloc() {
    genosig_t* sig;
    sig = (genosig_t*)calloc(1, sizeof(genosig_t));
    return sig;
}
void genosig_free(genosig_t* target) {
    if (target) {
        genosig_clearsig(target);
        genosig_cleargenomedata(target);
        if(target->allocatedname) free(target->name);
        if (target->lineage) free(target->lineage);
        memset(target, 0, sizeof(genosig_t));
        free(target);
    }
}

size_t genosig_info_sigsize(genosig_t* target) { return target->sigsize; }
size_t genosig_info_length(genosig_t* target) { return target->siglen; }
char* genosig_info_name(genosig_t* target) { return target->name; }
int genosig_info_ismuxed(genosig_t* target) { return (target->sigtype == _SIGTYPE_MUXED); }
int genosig_matchesname(genosig_t* target, int64_t taxid, char* literalname) {
    int result;
    result = 0;
    if (!target) return 0;
    if (target->lineage) {
        result = lineage_includes_id(target->lineage, taxid);
    }
    if (!result && target->filepath) {
        if (strcmp(literalname, target->filepath)==0)result = 1;
    }
    if (!result && target->name) {
        if (strcmp(literalname, target->name)==0)result = 1;
    }
    if (!result && target->filepath) {
        if (beginswith(literalname, target->filepath))result = 1;
    }
    if (!result && target->name) {
        if (beginswith(literalname, target->name))result = 1;
    }
    return result;
}
/* import/export data */
genosig_t* genosig_fullgenome(nucseq** data, size_t nseqs, int is_single_strand, int copylvl){
    genosig_t* res;
    nucseq** data2;
    size_t i;
    res = (genosig_t*)calloc(1,sizeof(genosig_t));
    if (copylvl == COPYLVL_NOTHING) {
        if (data && nseqs > 0)
            res->linkedgenome = 1;
    }
    else if (copylvl == COPYLVL_INTEGRATE) {
        if (data && nseqs > 0)
            res->linkedgenome = 0;
    }
    else if (copylvl == COPYLVL_FULL) {
        data2 = (nucseq**) calloc(nseqs, sizeof(nucseq*));
        for (i = 0;i < nseqs;i++) {
            data2[i] = (nucseq*) calloc(1, sizeof(nucseq));
            subseq(data2[i], data[i], 0, data[i]->len);
        }
        data = data2;
    }
    res->contigs = data;
    res->ncontigs = nseqs;
    res->signame = _SIGNAME_SEQUENCE;
    res->sigtype = _SIGTYPE_NOSIG;
    res->weight = 1.0;
    res->singlestrand = is_single_strand;
    return res;
}
genosig_t* genosig_linkfilename(genosig_t* sig, char* fn) {
    sig->filepath = fn;
    return sig;
}
genosig_t* genosig_genomefileref(char* path) {
    return genosig_linkfilename(genosig_alloc(), path);
}
genosig_t* genosig_linkname(genosig_t* sig, char* name) {
    if (sig->allocatedname) free(sig->name);
    sig->name = name;
    sig->allocatedname = 0;
    return sig;
}
genosig_t* genosig_autoname(genosig_t* sig, char* filepath) {
    char* basefn;
    size_t idend;
    if (sig->name && sig->allocatedname) free(sig->name);
    sig->filepath = filepath;
    basefn = os_rmdirname(filepath);
    if (beginswith("GCF_", basefn) || beginswith("GCA_", basefn)) {
        idend = text_findchar(basefn + 5, '_') + 5;
        sig->name = malloc(idend + 1);
        memcpy(sig->name, basefn, idend);
        sig->name[idend] = 0;
        sig->allocatedname = 1;
    }
    else {
        idend = strlen(basefn);
        sig->name = malloc(idend + 1);
        memcpy(sig->name, basefn, idend);
        sig->name[idend] = 0;
        sig->allocatedname = 1;
    }
    return sig;
}
genosig_t* genosig_renameby(genosig_t* sig, taxtree_t* tree, int copylvl, int mode) {
    int64_t id;
    size_t tmplen;
    char* tmpptr;
    id = -1;
    if (sig->filepath) id = taxtree_getid(tree, sig->filepath);
    if (id <= 0 && sig->name) id = taxtree_getid(tree, sig->name);
    if (id >= 0) {
        tmpptr = NULL;
        if (mode == RENAMEBY_NODENAME) tmpptr = taxtree_getnameptr(tree, id);
        if (mode == RENAMEBY_NODELABEL) tmpptr = taxtree_getlabelptr(tree, id);
        if (tmpptr) {
            if (sig->allocatedname) free(sig->name);
            if (copylvl == COPYLVL_NOTHING) {
                sig->name = tmpptr;
                sig->allocatedname = 0;
            }
            else if(copylvl == COPYLVL_INTEGRATE) {
                sig->name = tmpptr;
                sig->allocatedname = 1;
            }
            else if (copylvl == COPYLVL_FULL) {
                tmplen = strlen(tmpptr);
                sig->name = malloc(tmplen + 1);
                memcpy(sig->name, tmpptr, tmplen+1);
                sig->allocatedname = 1;
            }
        }
        if (mode == RENAMEBY_SPECIES) {
            tmpptr = taxtree_generatereadablename(tree, id);
            if (tmpptr) {
                sig->name = tmpptr;
                sig->allocatedname = 1;
            }
        }
    }
    return sig;
}
genosig_t* genosig_integratelineage(genosig_t* sig, lineage_t* lineage) {
    sig->lineage = lineage;
    return sig;
}
genosig_t* genosig_importlineage(genosig_t* sig, taxtree_t* tree, int* nullflag) {
    int64_t id;
    id = -1;
    if (sig->filepath) id = taxtree_getid(tree, sig->filepath);
    if (id <= 0 && sig->name) id = taxtree_getid(tree, sig->name);
    sig->lineage = taxtree_lineageof(tree, id);
    if (sig->lineage) *nullflag = 0;
    else *nullflag = 1;
    return sig;
}
genosig_t* genosig_makewindows(genosig_t* target, size_t winsize, size_t winskip){
    target->window = winsize;
    target->windowskip = winskip;
    return target;
}
genosig_t* genosig_suffixtrees(genosig_t* sig) {
    size_t seqid;
    suftree_t** res;
    res = (suftree_t**)malloc(sizeof(suftree_t*)*sig->ncontigs);
    for (seqid = 0;seqid < sig->ncontigs;seqid++) {
        res[seqid] = suftree_from(sig->contigs[seqid]->seq, SUFTREE_DEFAULT, sig->contigs[seqid]->len);
    }
    sig->sigdata = res;
    sig->siglen = sig->ncontigs;
    sig->sigsize = sizeof(suftree_t*)*sig->ncontigs;
    sig->sigtype = _SIGTYPE_SUFTREEARRAY;
    sig->signame = _SIGNAME_SUFTREE;
    return sig;
}
static inline suftree_t* _genosig_findgene(nucseq* data, suftree_t* gene, size_t minlen, size_t bridgegap, size_t minseed, double* alnscore) {
    char* marks;
    char** detected;
    char** detectedrc;
    double* alignment_score;
    size_t* len;
    size_t* rclen;
    size_t i, n, nrc, maxlenid, resultlen;
    suftree_t* result;
    nucseq tmpseq = EMPTYSEQ;

    resultlen = 0;
    alignment_score = NULL;
    detected = NULL;
    len = NULL;
    n = nrc = 0;
    *alnscore = -1.0;
    marks = suftree_approximatesearch(data->seq, data->len, gene, minseed, minlen, bridgegap, 1);
    detected = marks_split(data->seq, marks, data->len, &len, &n, KEEP_MARKED);
    nucseqrevcomp(data, &tmpseq);
    free(marks);
    marks = suftree_approximatesearch(tmpseq.seq, tmpseq.len, gene, minseed, minlen, bridgegap, 1);
    detectedrc = marks_split(tmpseq.seq, marks, tmpseq.len, &rclen, &nrc, KEEP_MARKED);
    free(marks);
    clear_nucseq(&tmpseq);
    if (n == 0 && nrc > 0) {
        if (detected)free(detectedrc);
        if (len) free(len);
        n = nrc;
        detected = detectedrc;
        len = rclen;
    }
    else if (n > 0 && nrc > 0) {
        detected = (char**)realloc(detected, sizeof(char*)*(n + nrc));
        len = (size_t*)realloc(len, sizeof(size_t)*(n + nrc));
        memcpy(detected + n, detectedrc, sizeof(char*)*nrc);
        memcpy(len + n, rclen, sizeof(size_t)*nrc);
        free(detectedrc);
        free(rclen);
        n += nrc;
    }
    if (n > 0) {
        alignment_score = (double*)malloc(sizeof(double)*n);
        alignment_score[0] = 1.0;
        maxlenid = 0;
        for (i = 0;i < n;i++) {
            alignment_score[i] = suftree_roughalign(detected[i], len[i], gene, 8, 1, 0, NULL);
            if (alignment_score[i] < alignment_score[maxlenid])maxlenid = i;
        }
        result = suftree_from(detected[maxlenid], SUFTREE_DEFAULT | SUFTREECOPYTEXT, len[maxlenid]);
        *alnscore = 1.0 - alignment_score[maxlenid];
        resultlen = len[maxlenid];
    }
    else {
        result = NULL;
        resultlen = 0;
    }
    for (i = 0;i < n;i++) {
        free(detected[i]);
    }
    if (alignment_score)
        free(alignment_score);
    if (detected)
        free(detected);
    if (len)
        free(len);
    return result;
}
genosig_t* genosig_findgenes(genosig_t* data, int is_single_strand, genosig_t* genesigs, int detectionmethod) {
    size_t geneid, seqid;
    suftree_t* tmptree;
    suftree_t** output;
    suftree_t** genes;
    double alnscore;
    double* alnscores;
    size_t bridgegap;
    size_t minlen;
    size_t minseed;
    
    if (genesigs->sigtype != _SIGTYPE_SUFTREEARRAY) return data;
    genes = (suftree_t**)(genesigs->sigdata);
    output = calloc(genesigs->siglen, sizeof(suftree_t*));
    alnscores = calloc(genesigs->siglen, sizeof(double));
    for (geneid = 0; geneid < genesigs->siglen; geneid++) {
        minlen = suftree_datalen(genes[geneid])*5/8;
        minseed = 9 + minlen / 300;
        bridgegap = minlen / 4;
        for (seqid = 0; seqid < data->ncontigs; seqid++) {
            tmptree = _genosig_findgene(data->contigs[seqid], genes[geneid], minlen, bridgegap, minseed, &alnscore);
            if (alnscore > alnscores[geneid]) {
                alnscores[geneid] = alnscore;
                if (output[geneid])suftree_free(output[geneid]);
                output[geneid] = tmptree;
            }
            else if(tmptree)
                suftree_free(tmptree);
        }
    }
    genosig_clearsig(data);
    data->sigdata = output;
    data->siglen = genesigs->siglen;
    data->sigsize = genesigs->siglen * sizeof(suftree_t*);
    data->sigtype = _SIGTYPE_SUFTREEARRAY;
    data->signame = _SIGNAME_SUFTREE;

    return data;
}
genosig_t* genosig_loadtxt(PF_t* file){
    return NULL;
}
genosig_t* genosig_loadbin(PF_t* file){
    genosig_t* out;
    genosig_t* tmp;
    char* tmpstr;
    uint16_t unused;
    size_t i;
    uint64_t labellen;
    int incorrect;
    out = genosig_alloc();
    unused = PFgetint16(file);
    out->signame = PFgetint16(file);
    out->sigtype = PFgetint16(file);
    out->siglen = PFgetint64(file);
    out->sigsize = PFgetint64(file);
    labellen = PFgetint64(file);
    out->name = malloc(labellen + 1);
    if (!(out->name) && labellen > 0) {
        genosig_free(out);
        return NULL;
    }
    i = PFread(out->name, 1, labellen, file);
    if (i != labellen) {
        genosig_free(out);
        return NULL;
    }
    out->name[i] = 0;
    out->allocatedname = 1;
    /* sanity check */
    if (out->signame >= _SIGNAME_MAX
        || out->sigtype >= _SIGTYPE_MAX
        || out->siglen<0
        || out->sigsize < out->siglen
        || strlen(out->name)<labellen) {
        genosig_free(out);
        return NULL;
    }
    incorrect = 0;
    switch (out->signame) {
    case _SIGNAME_SEQUENCE:
        out->ncontigs = out->siglen;
        out->siglen = 0;
        out->sigsize = 0;
        out->contigs = calloc(out->ncontigs, sizeof(nucseq*));
        for (i = 0;i < out->ncontigs;i++) {
            tmpstr = PFreadunimem64(file);
            if (!tmpstr) {
                incorrect = 1;
                out->ncontigs = i;
                break;
            }
            out->contigs[i] = calloc(1, sizeof(nucseq));
            nucseq_fromunimem(out->contigs[i],tmpstr);
            free(tmpstr);
        }
        break;
    case _SIGNAME_COUNTS:
    case _SIGNAME_FREQPROF:
    case _SIGNAME_MMZ:
    case _SIGNAME_KARLIN:
    case _SIGNAME_HASHLIST:
    case _SIGNAME_PRESABS:
    case _SIGNAME_GC:
    case _SIGNAME_GENOMELEN:
    case _SIGNAME_KARLINL:
    case _SIGNAME_MULTINORM:
        out->sigdata = malloc(out->sigsize);
        i = PFread(out->sigdata, 1, out->sigsize, file);
        if (i < out->sigsize) incorrect = 1;
        break;
    case _SIGNAME_SUFTREE:
        break;
    case _SIGNAME_MUXED:
        out->sigdata = (genosig_t**)malloc(out->sigsize * sizeof(genosig_t*));
        for (i = 0;i < out->siglen;i++) {
            tmp = genosig_loadbin(file);
            if (!tmp) {
                incorrect = 1;
                out->siglen = i;
                break;
            }
            ((genosig_t**)(out->sigdata))[i] = tmp;
        }
        break;
    default:
        incorrect = 1;
        break;
    }
    if (incorrect) {
        genosig_free(out);
        out = NULL;
    }
    return out;
}
char* _genosig_astxt(genosig_t* sig, size_t* outsize, char* prefixstr, char separator) {
    char* result;
    char* tmpptr;
    char** namecomponents;
    char** resultparts;
    char* signame;
    char* suffix;
    size_t* partlens;
    size_t nparts;
    size_t tmpsz, signamelen, siglen;
    size_t resultlen;
    size_t i;
    int64_t valuei;
    double valued;
    int issinglevec;
    nucseq tmpseq = EMPTYSEQ;
    signame = "";
    partlens = NULL;
    resultparts = NULL;
    nparts = 0;
    if (sig->signame == _SIGNAME_NONE) {
        *outsize = 0;
        return NULL;
    }
    switch (sig->signame) {
        case _SIGNAME_SEQUENCE:
            signame = "Sequence";
            break;
        case _SIGNAME_COUNTS:
            signame = "KmerCount";
            break;
        case _SIGNAME_FREQPROF:
            signame = "FrequencyProfile";
            break;
        case _SIGNAME_MMZ:
            signame = "MarkovModelZScore";
            break;
        case _SIGNAME_KARLIN:
        case _SIGNAME_KARLINL:
            signame = "AltKarlinSignature";
            break;
        case _SIGNAME_MULTINORM:
            signame = "NDimNormalParameters";
            break;
        case _SIGNAME_HASHLIST:
            signame = "ListOfHashes";
            break;
        case _SIGNAME_PRESABS:
            signame = "PresenceAbsence";
            break;
        case _SIGNAME_SUFTREE:
            signame = "Sequence";
            break;
        case _SIGNAME_MUXED:
            signame = "group";
            break;
        case _SIGNAME_FILENAME:
            signame = NULL;
        default:
            break;
    }
    issinglevec = 0;
    resultlen = 0;
    switch (sig->sigtype) {
        case _SIGTYPE_NOSIG:
            resultlen = 0;
            resultparts = (char**)malloc(sizeof(char*)*sig->ncontigs);
            partlens = (size_t*)malloc(sizeof(size_t)*sig->ncontigs);
            nparts = sig->ncontigs;
            for (i = 0;i < sig->ncontigs;i++) {
                resultparts[i] = nucseq2fasta(sig->contigs[i], NULL);
                tmpsz = strlen(resultparts[i]);
                partlens[i] = tmpsz;
                resultlen += tmpsz;
            }
            break;
        case _SIGTYPE_BYTE:
        case _SIGTYPE_I32:
        case _SIGTYPE_I64:
        case _SIGTYPE_DBL:
        case _SIGTYPE_SZT:
            issinglevec = 1;
            break;
        case _SIGTYPE_SUFTREEARRAY:
        case _SIGTYPE_HASHTBL:
        case _SIGTYPE_MUXED:
            break;
        default:
            break;
    }
    if (issinglevec) {
        /* all single vector signatures follow roughly the same pattern */
        siglen = sig->siglen;
        if (sig->signame == _SIGNAME_MULTINORM) {
            siglen *= 2;
        }
        resultparts = (char**)malloc(sizeof(char*)*(siglen + 1));
        partlens = (size_t*)malloc(sizeof(size_t)*(siglen + 1));
        nparts = siglen + 1;
        /* signature name */
        namecomponents = calloc(3,sizeof(char*));
        i = 0;
        if (prefixstr) {
            namecomponents[i] = prefixstr;
            i++;
        }
        if (sig->signame == _SIGNAME_FILENAME) {
            namecomponents[i] = sig->filepath;
            i++;
        }
        else if (sig->name) {
            namecomponents[i] = sig->name;
            i++;
        }
        if (signame) {
            namecomponents[i] = signame;
            i++;
        }
        resultparts[0] = text_join(namecomponents, ".", NULL, NULL, i);
        partlens[0] = strlen(resultparts[0]);
        free(namecomponents);
        /* actual signature */
        switch (sig->sigtype) {
            case _SIGTYPE_BYTE:
                for (i = 0;i < siglen;i++) {
                    valuei = (int64_t)(((char*)(sig->sigdata))[i]);
                    resultparts[i + 1] = text_fromint(valuei);
                    partlens[i + 1] = strlen(resultparts[i + 1]);
                }
                break;
            case _SIGTYPE_I32:
                for (i = 0;i < siglen;i++) {
                    valuei = (int64_t)(((int32_t*)(sig->sigdata))[i]);
                    resultparts[i + 1] = text_fromint(valuei);
                    partlens[i + 1] = strlen(resultparts[i + 1]);
                }
                break;
            case _SIGTYPE_I64:
                for (i = 0;i < siglen;i++) {
                    valuei = ((int64_t*)(sig->sigdata))[i];
                    resultparts[i + 1] = text_fromint(valuei);
                    partlens[i + 1] = strlen(resultparts[i + 1]);
                }
                break;
            case _SIGTYPE_DBL:
                for (i = 0;i < siglen;i++) {
                    valued = ((double*)(sig->sigdata))[i];
                    resultparts[i + 1] = text_fromdbl(valued, 8);
                    partlens[i + 1] = strlen(resultparts[i + 1]);
                }
                break;
            case _SIGTYPE_SZT:
                for (i = 0;i < siglen;i++) {
                    valuei = (int64_t)(((size_t*)(sig->sigdata))[i]);
                    resultparts[i + 1] = text_fromint(valuei);
                    partlens[i + 1] = strlen(resultparts[i + 1]);
                }
                break;
            default:
                break;
        }
        
        for (i = 0;i < nparts;i++) {
            resultlen += partlens[i];
        }
        resultparts[nparts-1] = realloc(resultparts[nparts-1], partlens[nparts-1]+2);
        resultparts[nparts-1][partlens[nparts-1]] = '\n';
        partlens[nparts-1]++;
        resultlen++;
        resultparts[nparts - 1][partlens[nparts - 1]] = 0;
        resultlen += sizeof(char)*siglen; /* separators */
    }
    else if (sig->sigtype == _SIGTYPE_SUFTREEARRAY) {
        /* since suffix trees are somewhat difficult to visualise, output the sequences */
        resultparts = (char**)malloc(sizeof(char*)*(sig->siglen));
        partlens = (size_t*)malloc(sizeof(size_t)*(sig->siglen));
        nparts = sig->siglen;
        /* signature name */
        namecomponents = calloc(3, sizeof(char*));
        i = 0;
        if (prefixstr) {
            namecomponents[i] = prefixstr;
            i++;
        }
        if (sig->name) {
            namecomponents[i] = sig->name;
            i++;
        }
        namecomponents[i] = signame;
        i++;
        signame = text_join(namecomponents, ".", NULL, NULL, i);
        signamelen = strlen(signame);
        free(namecomponents);
        for (i = 0;i < sig->siglen;i++) {
            /* combine different things to create a nucseq */
            tmpseq.seq = suftree_getstrptr(((suftree_t**)(sig->sigdata))[i], &tmpsz);
            tmpseq.len = tmpsz;
            suffix = text_fromint((int64_t)i);
            tmpsz = strlen(suffix);
            tmpseq.name = malloc(signamelen + tmpsz + 2);
            memcpy(tmpseq.name, signame, signamelen);
            tmpseq.name[signamelen] = '.';
            memcpy(tmpseq.name + signamelen + 1, suffix, tmpsz);
            tmpseq.name[signamelen + tmpsz + 1] = 0;
            resultparts[i] = nucseq2fasta(&tmpseq, tmpseq.name);
            partlens[i] = strlen(resultparts[i]);
            /* free the name since we don't need it anymore */
            free(tmpseq.name);
        }
        for (i = 0;i < sig->siglen;i++) {
            resultlen += partlens[i];
        }
        free(signame);
    }
    result = (char*)malloc(resultlen+1);
    tmpptr = result;
    
    for (i = 0;i < nparts;i++) {
        memcpy(tmpptr, resultparts[i], partlens[i]);
        free(resultparts[i]);
        tmpptr += partlens[i];
        if (issinglevec && i < nparts-1) {
            *tmpptr = separator;
            tmpptr++;
        }
    }
    free(resultparts);
    free(partlens);
    result[resultlen] = 0;
    *outsize = resultlen;
    return result;
}
char* genosig_astxt(genosig_t* sig, size_t* outsize){
    return _genosig_astxt(sig, outsize, NULL, ',');
}
char* genosig_asbin(genosig_t* sig, size_t* outsize){
    char* result;
    char* tmpptr;
    char* label;
    char** resultparts;
    size_t* partlens;
    size_t nparts;
    size_t tmpsz;
    size_t resultlen;
    size_t i;
    uint64_t siglabellen;
    if (sig->signame == _SIGNAME_NONE) {
        *outsize = 0;
        return NULL;
    }
    resultparts = NULL;
    partlens = NULL;
    switch (sig->signame) {
        case _SIGNAME_SEQUENCE:
            resultlen = 0;
            resultparts = (char**) malloc(sizeof(char*)*sig->ncontigs);
            partlens = (size_t*)malloc(sizeof(size_t)*sig->ncontigs);
            nparts = sig->ncontigs;
            for (i = 0;i < sig->ncontigs;i++) {
                resultparts[i] = nucseq_unimem(sig->contigs[i], &tmpsz);
                partlens[i] = tmpsz;
                resultlen += tmpsz;
            }
            break;
        case _SIGNAME_COUNTS:
        case _SIGNAME_FREQPROF:
        case _SIGNAME_MMZ:
        case _SIGNAME_KARLIN:
        case _SIGNAME_HASHLIST:
        case _SIGNAME_PRESABS:
        case _SIGNAME_GC:
        case _SIGNAME_GENOMELEN:
        case _SIGNAME_KARLINL:
        case _SIGNAME_MULTINORM:
            resultlen = sig->sigsize;
            resultparts = (char**) malloc(sizeof(char*));
            partlens = (size_t*)malloc(sizeof(size_t));
            resultparts[0] = malloc(resultlen);
            memcpy(resultparts[0], sig->sigdata, sig->sigsize);
            partlens[0] = sig->sigsize;
            nparts = 1;
            break;
        case _SIGNAME_SUFTREE:
            resultlen = 0;
            nparts = 0;
            break;
        case _SIGNAME_MUXED:
            resultlen = 0;
            resultparts = (char**)malloc(sizeof(char*)*sig->siglen);
            partlens = (size_t*)malloc(sizeof(size_t)*sig->siglen);
            nparts = sig->siglen;
            for (i = 0;i < sig->siglen;i++) {
                resultparts[i] = genosig_asbin(((genosig_t**)(sig->sigdata))[i], &tmpsz);
                partlens[i] = tmpsz;
                resultlen += tmpsz;
            }
            break;
        default:
            resultlen = 0;
            nparts = 0;
            break;
    }
    resultlen += sizeof(uint16_t); /* placeholder */
    resultlen += sizeof(uint16_t); /* signature name */
    resultlen += sizeof(uint16_t); /* signature type */
    resultlen += sizeof(uint64_t); /* signature length */
    resultlen += sizeof(uint64_t); /* signature size */
    resultlen += sizeof(uint64_t); /* length of signature label */
    if (sig->name) label = sig->name;
    else if (sig->filepath) label = sig->filepath;
    else label = NULL;
    if (label)
        siglabellen = strlen(label);
    else
        siglabellen = 0;
    resultlen += sizeof(char)*siglabellen; /* signature label */

    result = (char*)malloc(resultlen);
    tmpptr = result;
    *((uint16_t*)tmpptr) = 0;
    tmpptr += sizeof(uint16_t);
    *((uint16_t*)tmpptr) = sig->signame;
    tmpptr += sizeof(uint16_t);
    *((uint16_t*)tmpptr) = sig->sigtype;
    tmpptr += sizeof(uint16_t);
    *((uint64_t*)tmpptr) = sig->siglen;
    tmpptr += sizeof(uint64_t);
    *((uint64_t*)tmpptr) = sig->sigsize;
    tmpptr += sizeof(uint64_t);
    *((uint64_t*)tmpptr) = siglabellen;
    tmpptr += sizeof(uint64_t);
    memcpy(tmpptr, label, siglabellen);
    tmpptr += sizeof(char)*siglabellen;
    for (i = 0;i < nparts;i++) {
        memcpy(tmpptr, resultparts[i], partlens[i]);
        free(resultparts[i]);
        tmpptr += partlens[i];
    }
    free(resultparts);
    free(partlens);
    *outsize = resultlen;
    return result;
}
void genosig_savetxt(genosig_t* sig, PF_t* file) {
    char* towrite;
    size_t fsize;
    towrite = genosig_astxt(sig, &fsize);
    if (towrite) {
        PFwrite(towrite, 1, fsize, file);
        free(towrite);
    }
}
void genosig_savebin(genosig_t* sig, PF_t* file) {
    char* towrite;
    size_t fsize;
    towrite = genosig_asbin(sig, &fsize);
    if (towrite) {
        PFwrite(towrite, 1, fsize, file);
        free(towrite);
    }
}
void genosig_export(genosig_t* sig, const char* filename) {
    PF_t* f;
    uint32_t prefix;
    PFopen(&f, filename, "wb");
    if (f) {
        prefix = (((uint32_t)'G') << 24) + (((uint32_t)'S') << 16);
        PFputint32(f, prefix);
        genosig_savebin(sig, f);
        PFclose(f);
    }
}
genosig_t* genosig_import_v1(PF_t* file) {
    /* old file format - type == 1 */
    genosig_t* result;
    int32_t datatype;
    uint64_t n;
    int64_t fnamelen;

#ifdef DEF_ARGPARSER
    args_report_warning(NULL, "File %s is stored in an old format - consider updating it.\n", PFgetbasenameptr(file));
#else
    fprintf(stderr, "File %s is stored in an old format - consider updating it.\n", PFgetbasenameptr(file));
#endif

    result = genosig_alloc();
    /* it is assumed that the file's endianness matches that of the OS */
    PFread(&datatype, sizeof(int32_t), 1, file);
    PFread(&n, sizeof(uint64_t), 1, file);
    if (datatype == 2) {
        /* rem: datatype == 2 is bytes */
        result->sigdata = malloc(n);
        result->sigtype = _SIGTYPE_BYTE;
        result->signame = _SIGNAME_UNKNOWN;
        PFread(result->sigdata, 1, n, file);
    }
    else {
        result->sigdata = malloc(sizeof(double)*n);
        result->sigtype = _SIGTYPE_DBL;
        result->signame = _SIGNAME_UNKNOWN;
        PFread(result->sigdata, sizeof(double), n, file);
    }
    result->siglen = n;
    fnamelen = PFgetint64(file);
    if (fnamelen > 0) {
        result->name = malloc(sizeof(char)*(fnamelen + 1));
        result->name[fnamelen] = 0;
        result->allocatedname = 1;
        PFread(result->name, sizeof(char), fnamelen, file);
    }
    return result;
}
genosig_t* genosig_import_v1muxed(PF_t* file) {
    genosig_t** allsigs;
    genosig_t* output;
    int64_t genocount;
    int32_t irrelevant;
    int64_t i;
    genocount = PFgetint64(file);
    allsigs = (genosig_t**)calloc(genocount, sizeof(genosig_t*));
    if (!allsigs) return NULL;
    for (i = 0;i < genocount;i++) {
        irrelevant = PFgetint32(file);
        allsigs[i] = genosig_import_v1(file);
    }
    output = genosig_multiplex(allsigs, genocount);
    return output;
}

genosig_t* genosig_importextra(char* filename, uint32_t flags) {
    /*
    genosig_import will always consider fasta-formatted data as double-stranded.
    Use genosig_fullgenome(..., ..., 1, ...) when dealing with single-stranded
    data instead of genosig_import.
    */
    int32_t formattype;
    genosig_t* result;
    genosig_t** tmp;
    nucseq** nsarr;
    size_t i, count;
    PF_t* f;
    int32_t prefix1;
    int32_t prefix2;
    size_t badcount;
    prefix1 = (((uint32_t)'G') << 24) + (((uint32_t)'S') << 16);
    prefix2 = (((uint32_t)'S') << 8) + ((uint32_t)'G');
    result = NULL;
    PFopen(&f, filename, "rb");
    if (f) {
        if (!endswith(".sig", filename) && !endswith(".sdb", filename)) {
            nsarr = nucseq_array_from_fasta(f, &count, 1, 2, &badcount);
            result = genosig_fullgenome(nsarr, count, 0, COPYLVL_INTEGRATE);
        }
        else {
            formattype = PFgetint32(f);
            if (formattype == 1) {
                result = genosig_import_v1(f);
            }
            else if (PFseemslikeASCIItext(f, 200)) {
                if (0) /* (flags&GENOSIG_IMPORT_NOFASTA) */ {
                    result = NULL;
                }
                else {
                    PFrewind(f);
                    nsarr = nucseq_array_from_fasta(f, &count, 1, 1, &badcount);
                    result = genosig_fullgenome(nsarr, count, 0, COPYLVL_INTEGRATE);
                    tmp = genosig_demux(result, &count);
                    for (i = 0;i < count;i++) {
                        tmp[i] = genosig_suffixtrees(tmp[i]);
                        genosig_cleargenomedata(tmp[i]);
                    }
                    result = genosig_multiplex(tmp, count);
                }
            }
            else if (formattype == prefix1 || formattype == prefix2) {
                result = genosig_loadbin(f);
            }
            if (!result && endswith(".sdb", filename)) {
                PFrewind(f);
                result = genosig_import_v1muxed(f);
            }
        }
        if (result && !(result->name)) {
            genosig_autoname(result, filename);
        }
        if (result && !(result->filepath))
            result->filepath = os_rmdirname(filename);
        PFclose(f);
    }
    return result;
}
genosig_t* genosig_import(char* filename) {
    return genosig_importextra(filename, 0);
}
/*
about multiplexed signatures
any signature involved in a (de)multiplexing operation should only be free'd once.
As such, signatures should not free'd after multiplexing (i.e. should not be free'd individually).
On the other hand, demultiplexing a signature invalidates its pointer.
*/
genosig_t* genosig_multiplex(genosig_t** sigs, size_t count){
    genosig_t* result;
    result = genosig_alloc();
    result->sigtype = _SIGTYPE_MUXED;
    result->signame = _SIGNAME_MUXED;
    result->siglen = count;
    result->sigsize = sizeof(genosig_t*)*count;
    result->sigdata = sigs;
    result->sigmux = 1;
    return result;
}
genosig_t** genosig_demux(genosig_t* muxedsigs, size_t* count){
    genosig_t** output;
    nucseq** tmp;
    nucseq** dataptr;
    size_t i,lcount;
    output = NULL;
    if (muxedsigs->sigtype == _SIGTYPE_MUXED) {
        lcount = muxedsigs->siglen;
        genosig_cleargenomedata(muxedsigs);
        muxedsigs->sigsize = 0;
        muxedsigs->sigtype = _SIGTYPE_NOSIG;
        muxedsigs->signame = _SIGNAME_NONE;
        muxedsigs->siglen = 0;
        output = (genosig_t**)muxedsigs->sigdata;
        free(muxedsigs);
    }
    else if (muxedsigs->window > 100 && muxedsigs->windowskip > 0) {
        tmp = nucseq_winreadmultiple(muxedsigs->contigs, muxedsigs->ncontigs, &lcount, muxedsigs->window, muxedsigs->windowskip);
        output = (genosig_t**)malloc(sizeof(genosig_t*)*lcount);
        for (i = 0;i < lcount;i++) {
            dataptr = (nucseq**)malloc(sizeof(nucseq*));
            dataptr[0] = tmp[i];
            output[i] = genosig_fullgenome(dataptr, 1, 0, COPYLVL_INTEGRATE);
            output[i]->weight = (double)(dataptr[0]->len);
        }
        free(tmp);
        genosig_free(muxedsigs);
    }
    else if (muxedsigs->ncontigs > 1) {
        tmp = muxedsigs->contigs;
        lcount = muxedsigs->ncontigs;
        output = (genosig_t**)malloc(sizeof(genosig_t*)*lcount);
        for (i = 0;i < lcount;i++) {
            if (muxedsigs->linkedgenome) {
                output[i] = genosig_fullgenome(tmp + i, 1, 0, COPYLVL_NOTHING);
                output[i]->linkedgenome = muxedsigs->linkedgenome;
                if (tmp[i]->name) {
                    output[i]->name = memcpyalloc(tmp[i]->name, strlen(tmp[i]->name)+1);
                    output[i]->allocatedname = 1;
                }
                output[i]->weight = (double)(tmp[i]->len);
            }
            else {
                dataptr = (nucseq**)malloc(sizeof(nucseq*));
                dataptr[0] = tmp[i];
                output[i] = genosig_fullgenome(dataptr, 1, 0, COPYLVL_NOTHING);
                output[i]->linkedgenome = 0;
                if (tmp[i]->name) {
                    output[i]->name = memcpyalloc(tmp[i]->name, strlen(tmp[i]->name)+1);
                    output[i]->allocatedname = 1;
                }
                output[i]->weight = (double)(dataptr[0]->len);
            }
        }
        if (!(muxedsigs->linkedgenome)) {
            free(tmp);
            muxedsigs->contigs = NULL;
        }
        muxedsigs->linkedgenome = 1;
        genosig_free(muxedsigs);
    }
    *count = lcount;
    return output;
}

/* basis */
genosig_t* genosig_keepfilenameonly(genosig_t* sig, size_t unused) {
    genosig_cleargenomedata(sig);
    sig->signame = _SIGNAME_FILENAME;
    sig->sigtype = _SIGTYPE_BYTE;
    genosig_autoname(sig, sig->filepath);
    return sig;
}
genosig_t* genosig_GC(genosig_t* sig, size_t unused) {
    int64_t GC,ACGT;
    double result;
    size_t i,j;
    genosig_clearsig(sig);
    sig->sigsize = sizeof(double);
    sig->siglen = 1;
    GC = ACGT = 0;
    for (i = 0;i < sig->ncontigs;i++) {
        for (j = 0;j < sig->contigs[i]->len;j++) {
            if (sig->contigs[i]->seq[j] == nucC || sig->contigs[i]->seq[j] == nucG) GC++;
        }
        ACGT += sig->contigs[i]->len;
    }
    result = ((double)GC) / (double)ACGT;
    sig->sigdata = malloc(sig->sigsize);
    sig->sigtype = _SIGTYPE_DBL;
    sig->signame = _SIGNAME_GC;
    memcpy(sig->sigdata, &result, sig->sigsize);
    return sig;
}
genosig_t* genosig_length(genosig_t* sig, size_t unused){
    size_t i;
    size_t tmp;
    genosig_clearsig(sig);
    sig->sigsize = sizeof(size_t);
    sig->siglen = 1;
    tmp = 0;
    for (i = 0;i < sig->ncontigs;i++) {
        tmp += sig->contigs[i]->len;
    }
    sig->sigdata = malloc(sig->sigsize);
    sig->sigtype = _SIGTYPE_SZT;
    sig->signame = _SIGNAME_GENOMELEN;
    memcpy(sig->sigdata, &tmp, sig->sigsize);
    return sig;
}
static uint32_t* _count_short_kmers(nucseq* sequence, size_t k, size_t* outcount, int addrccounts) {
    int64_t curn = 0;
    size_t numvals = 1;
    size_t contiglen = 0;
    size_t skipamnt;
    int64_t kmer, kmerrc, firstnuc, curnuc, neg, neg0;
    size_t mask;
    uint32_t* counts = NULL;
    uint32_t* count2add = NULL;
    /* high lengths do not make sense with this implementation, but we assume that proper check have been done */
    /* count the number of cells we need to have. No point in using pow here. */
    numvals <<= (2LL * k);
    mask = numvals - 1;
    if (outcount) *outcount = numvals;
    counts = (uint32_t*)malloc(sizeof(uint32_t)*numvals);
    if (!counts) {
        fprintf(stderr, "Not enough memory for count");
    }
    memset(counts, 0, sizeof(uint32_t)*numvals);
    kmer = 0;
    skipamnt = k - 1;
    for (curn = 0;curn < (int64_t)(sequence->len);curn++) {
        if (sequence->seq[curn] >= 0)
        {
            contiglen++;
            kmer = ((kmer << 2LL) & mask);
            kmer += sequence->seq[curn];
            if (skipamnt)
                skipamnt--;
            else
                counts[kmer]++;
        }
        else
        {
            skipamnt = k - 1;
            contiglen = 0;
        }
    }
    if (addrccounts) {
        count2add = (uint32_t*)malloc(sizeof(uint32_t)*numvals);
        kmerrc = mask;
        neg0 = (numvals >> 2);
        firstnuc = numvals - (numvals >> 2);
        curnuc = firstnuc;
        for (kmer = 0;kmer < (int64_t)numvals;kmer++) {
            count2add[kmerrc] = counts[kmer];
            neg = neg0;
            curnuc = firstnuc;
            while (curnuc && (curnuc&kmerrc) == 0) {
                kmerrc += curnuc;
                curnuc >>= 2;
                neg >>= 2;
            }
            kmerrc -= neg;
        }
        veci32_add(counts, count2add, numvals);
        free(count2add);
    }
    return counts;
}
static void* _count_long_kmers(nucseq* sequence, size_t k) {
    fprintf(stderr, "Long k-mer counting not implemented yet");
    return NULL;
}
genosig_t* genosig_kmercount(genosig_t* sig, size_t k){
    uint32_t* tmpcount;
    size_t i;
    size_t ncomb;
    /* kmercount overwrites other signatures */
    genosig_clearsig(sig);
    if (k>0 && k < 16) {
        if (sig->ncontigs) {
            sig->sigdata = _count_short_kmers(sig->contigs[0], k, &ncomb, !(sig->singlestrand));
            sig->siglen = ncomb;
            sig->sigsize = ncomb * sizeof(uint32_t);
            sig->sigtype = _SIGTYPE_I32;
            sig->signame = _SIGNAME_COUNTS;
            for (i = 1;i < sig->ncontigs;i++) {
                tmpcount = (uint32_t*)_count_short_kmers(sig->contigs[i], k, NULL, !(sig->singlestrand));
                veci32_add(sig->sigdata, tmpcount, ncomb);
                free(tmpcount);
            }
            /*note: [contigs] contains linked data, and should therefore be free'd separately*/
        }
    }
    return sig;
}
#define MINMINHASHLEN   12
genosig_t* genosig_minhash(genosig_t* sig, size_t numhashes){
    int k;
    genosig_clearsig(sig);
    if (sig->window <= MINMINHASHLEN)k = 21;
    else k = (int)(sig->window);
    sig->sigdata = minhash_Msig(sig->contigs, sig->ncontigs, k, numhashes, 200000);
    sig->sigtype = _SIGTYPE_I64;
    sig->siglen = numhashes;
    sig->signame = _SIGNAME_HASHLIST;
    sig->sigsize = numhashes * sizeof(int64_t);
    return sig;
}
genosig_t* genosig_kmerfreq(genosig_t* sig, size_t k){
    double* freqs;
    uint64_t totcount;
    sig = genosig_kmercount(sig, k);
    totcount = veci32_sum64((int32_t*)(sig->sigdata), sig->siglen);
    freqs = veci32_to_vec((int32_t*)(sig->sigdata), sig->siglen);
    if (totcount > 0) {
        vec_scale(freqs, 1.0 / ((double)totcount), sig->siglen);
    }
    else {
        vec_scale(freqs, 0.0, sig->siglen);
    }
    free(sig->sigdata);
    sig->sigdata = freqs;
    sig->sigsize = sig->siglen * sizeof(double);
    sig->sigtype = _SIGTYPE_DBL;
    sig->signame = _SIGNAME_FREQPROF;
    return sig;
}
double* _markov_model_diff_zscore_infinite(int32_t* full, int32_t* leftright, int32_t* centre, size_t k) {
    double* freqs;
    double E;
    double var;
    double l, r, m;
    uint64_t rmod;
    int64_t N;
    uint64_t i;
    size_t siglen;
    freqs = NULL;
    if (k > 2) {
        siglen = 1LL << (2 * k);
        rmod = siglen >> 2;
        N = veci32_sum64(full, siglen);

        freqs = veci32_to_vec(full, siglen);
        for (i = 0; i < siglen; i++) {
            l = (double)(leftright[i >> 2]);
            r = (double)(leftright[i%rmod]);
            m = (double)(centre[(i%rmod) >> 2]);
            E = l * r / m; /* calculate expected value based on (k-1)-mers */
            var = E*(m - l)*(m - r) / (m*m); /* derived from the variance formula for bionomial distributions, I think */
            freqs[i] = (freqs[i] - E) / (sqrt(var)); /* z-score = (Value-expected)/stddev */
        }
    }
    return freqs;
}
genosig_t* genosig_mmz_inf(genosig_t* sig, size_t k) {
    int32_t* counts;
    int32_t* counts_1;
    int32_t* counts_2;
    int32_t* tmpcount;
    double* freqs;
    size_t ncomb;
    size_t i;
    if (k>2 && k < 16) {
        if (sig->ncontigs) {
            counts = _count_short_kmers(sig->contigs[0], k, &ncomb, !(sig->singlestrand));
            sig->siglen = ncomb;
            for (i = 1;i < sig->ncontigs;i++) {
                tmpcount = (uint32_t*)_count_short_kmers(sig->contigs[i], k, NULL, !(sig->singlestrand));
                veci32_add(counts, tmpcount, ncomb);
                free(tmpcount);
            }
            counts_1 = _count_short_kmers(sig->contigs[0], k-1, &ncomb, !(sig->singlestrand));
            for (i = 1;i < sig->ncontigs;i++) {
                tmpcount = (uint32_t*)_count_short_kmers(sig->contigs[i], k-1, NULL, !(sig->singlestrand));
                veci32_add(counts_1, tmpcount, ncomb);
                free(tmpcount);
            }
            counts_2 = _count_short_kmers(sig->contigs[0], k-2, &ncomb, !(sig->singlestrand));
            for (i = 1;i < sig->ncontigs;i++) {
                tmpcount = (uint32_t*)_count_short_kmers(sig->contigs[i], k-2, NULL, !(sig->singlestrand));
                veci32_add(counts_2, tmpcount, ncomb);
                free(tmpcount);
            }
            freqs = _markov_model_diff_zscore_infinite(counts, counts_1, counts_2, k);
            free(counts);
            free(counts_1);
            free(counts_2);
            genosig_clearsig(sig);
            sig->sigdata = freqs;
            sig->sigsize = sizeof(double)*sig->siglen;
            sig->sigtype = _SIGTYPE_DBL;
            sig->signame = _SIGNAME_MMZ;
        }
    }
    return sig;
}
double* _markov_model_diff_zscore(int32_t* counts, size_t k, size_t siglen) {
    double* freqs;
    double E;
    double var;
    double l, r, m;
    int32_t* lcount;
    int32_t* rcount;
    int32_t* mcount;
    uint64_t rmod;
    uint64_t i;
    freqs = NULL;
    if (k > 2) {
        lcount = calloc((siglen >> 2), sizeof(double));
        rcount = calloc((siglen >> 2), sizeof(double));
        mcount = calloc((siglen >> 4), sizeof(double));

        rmod = siglen >> 2;

        for (i = 0; i < siglen; i++) {
            lcount[i >> 2] += counts[i];
            rcount[i%rmod] += counts[i];
            mcount[(i%rmod) >> 2] += counts[i];
        }

        freqs = veci32_to_vec(counts, siglen);
        for (i = 0; i < siglen; i++) {
            l = (double)(lcount[i >> 2]);
            r = (double)(rcount[i%rmod]);
            m = (double)(mcount[(i%rmod) >> 2]);
            E = l * r / m; /* calculate expected value based on (k-1)-mers */
            var = E*(m - l)*(m - r) / (m*m); /* derived from the variance formula for bionomial distributions, I think */
            freqs[i] = (freqs[i] - E) / (sqrt(var)); /* z-score = (Value-expected)/stddev */
        }
        free(lcount);
        free(rcount);
        free(mcount);
    }
    return freqs;
}
genosig_t* genosig_mmz(genosig_t* sig, size_t k){
    double* freqs;
    sig = genosig_kmercount(sig, k);
    freqs = _markov_model_diff_zscore(sig->sigdata, k, sig->siglen);
    free(sig->sigdata);
    sig->sigdata = freqs;
    sig->sigsize = sizeof(double)*sig->siglen;
    sig->sigtype = _SIGTYPE_DBL;
    sig->signame = _SIGNAME_MMZ;
    return sig;
}
static inline int _count_set_bits(uint64_t value) {
/*use builtins when available */
#if GCC_VERSION > 30400
    return (int)__builtin_popcountll(value);
#elif (defined(_MSC_VER) && defined(_WIN64))
    return (int)__popcnt64(value);
#else
    /* For other ways to count bits, refer to:
    http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable  */
    uint64_t v = value;
    int c;
    for (c = 0; v; c++) {
        v &= v - 1; /* clear the least significant bit set */
    }
    return c;
#endif
}
static uint64_t _bindex(uint64_t S, uint64_t B, int setbits) {
    /* S represents [X_1 X_2 ... X_k] as SUM_n(X_n*4^n)
    /* B represents [B_1 B_2 ... B_k] as SUM_n(B_n*2^n)
    /* setbits is equal to SUM_n(B_n) */
    int c, X;
    uint64_t mask = 1;
    uint64_t div = 1;
    uint64_t mult = 1;
    uint64_t toret = -1;
    if (B != 0) {
        toret = 0;
        for (c = 0;c < setbits;c++) {
            /* count up to the next next non-zero bit */
            while ((mask&B) == 0) {
                mask <<= 1;
                div *= 4;
            }
            /* get the appropriate digit */
            X = (S / div) % 4;
            toret += X*mult;
            /* go to the next digit */
            mult *= 4;
            div *= 4;
            mask <<= 1;
        }
    }
    return toret;
}
static double* _karlinsig(uint32_t* counts, int k) {
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
    /* first calculate the number of elements in counts */
    /* number of S combinations: 4^k */
    if (!counts)return NULL;
    if (k < 1)return NULL;
    s_max = (uint64_t)pow(4, k);
    /* number of B combinations: 2^k-2 ( [0 0 ... 0] is not valid and [1 1 ... 1] is already stored in counts) */
    b_max = ((uint64_t)1 << k) - 2;
    /* allocate rho */
    rho = (double*)malloc(sizeof(double)*s_max);
    /* allocate CkSB */
    CkSB = (uint64_t**)malloc(sizeof(uint64_t*)*b_max);
    nCkSB = (int*)malloc(sizeof(int)*b_max);
    for (b = 0;b < b_max;b++) {
        nCkSB[b] = _count_set_bits(b + 1); /* 4^(count_set_bits(b+1)) */
        i = 1ULL << (2 * nCkSB[b]);
        CkSB[b] = (uint64_t*)malloc(sizeof(uint64_t)*i);
        memset(CkSB[b], 0, sizeof(uint64_t)*i);
    }

    /* Calculate CkSB (complexity: 8^k - original formula is more efficient for k>6 )*/
    /* iterate over all elements */
    for (s = 0; s < s_max; s++) {
        /* for each combination add the value to the proper place */
        for (b = 0;b < b_max;b++) {
            idx = _bindex(s, b + 1, nCkSB[b]);
            CkSB[b][idx] += counts[s];
        }
        CkS0 += counts[s];
    }
    /* calculate rho */
    for (s = 0;s < s_max;s++) {
        if (counts[s] > 0)
            rho[s] = (double)counts[s] / (double)CkS0;
        else
            rho[s] = 1.0 / (double)CkS0 / (double)s_max;
        for (b = 0;b < b_max;b++) {
            idx = _bindex(s, b + 1, nCkSB[b]);
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
    /* free CkSB and nCkSB - we do not need them anymore */
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
genosig_t* genosig_karlinsig(genosig_t* sig, size_t k){
    double* ksig;
    sig = genosig_kmercount(sig, k);
    ksig = _karlinsig((uint32_t*)(sig->sigdata), (int)k);
    if (sig->sigdata)
        free(sig->sigdata);
    if (ksig) {
        sig->sigdata = ksig;
        sig->sigsize = sig->siglen*(sizeof(double));
        sig->sigtype = _SIGTYPE_DBL;
        sig->signame = _SIGNAME_KARLIN;
    }
    else {
        sig->sigdata = ksig;
        sig->sigsize = 0;
        sig->sigtype = _SIGTYPE_NOSIG;
        sig->signame = _SIGNAME_NONE;
    }
    return sig;
}
genosig_t* genosig_karlinsigL(genosig_t* sig, size_t k) {
    double* ksig;
    double genolen;
    sig = genosig_kmercount(sig, k);
    genolen = (double)veci32_sum64((uint32_t*)(sig->sigdata), sig->siglen);
    ksig = _karlinsig((uint32_t*)(sig->sigdata), (int)k);
    free(sig->sigdata);
    if (ksig) {
        ksig = (double*)realloc(ksig, sizeof(double)*(sig->siglen + 1));
        ksig[sig->siglen] = genolen;
        sig->siglen++;
        sig->sigdata = ksig;
        sig->sigsize = sig->siglen*(sizeof(double));
        sig->sigtype = _SIGTYPE_DBL;
        sig->signame = _SIGNAME_KARLINL;
    }
    else {
        sig->sigtype = _SIGTYPE_NOSIG;
        sig->signame = _SIGNAME_NONE;
    }
    return sig;
}
genosig_t* genosig_presenceabsence(genosig_t* sig, size_t unused) {
    size_t i;
    size_t eltsize;
    size_t siglen;
    char* result;
    char* nullelt;
    char* tocompare;
    siglen = sig->siglen;
    eltsize = sig->sigsize / siglen;
    result = (char*) calloc(siglen, sizeof(char));
    nullelt = (char*)calloc(eltsize, sizeof(char));
    for (i = 0;i < sig->siglen;i++) {
        tocompare = (char*)(sig->sigdata) + i*eltsize;
        result[i] = (memcmp(tocompare, nullelt, eltsize) != 0);
    }
    /* note: the nature of the signature is not assumed, so it needs to be cleaned thouroughly */
    genosig_clearsig(sig);
    sig->sigdata = result;
    sig->siglen = siglen;
    sig->sigtype = _SIGTYPE_BYTE;
    sig->signame = _SIGNAME_PRESABS;
    return sig;
}
genosig_t* genosig_bactspeciessig(genosig_t* sig, size_t unused){
    return sig;
}
genosig_t* genosig_bacttypesig(genosig_t* sig, size_t unused){
    return sig;
}

genosig_t* genosig_avgsig(genosig_t* sig, size_t is_weighted) {
    /* average signature of a non-multiplexed signature is always itself */
    /* a new signature is generated otehrwise*/
    int16_t shared_sigtype;
    int16_t shared_signame;
    size_t shared_siglen;
    size_t i;
    genosig_t** sigarr;
    genosig_t* output;
    double* tmpd;
    double* tmpd2;
    double totweight;
    char* newname[2];
    if (sig->sigtype != _SIGTYPE_MUXED) return sig;
    sigarr = sig->sigdata;
    shared_sigtype = sigarr[0]->sigtype;
    shared_siglen = sigarr[0]->siglen;
    shared_signame = sigarr[0]->signame;
    for (i = 0;i < sig->siglen;i++) {
        if (sigarr[i]->sigtype != shared_sigtype)break;
        if (sigarr[i]->siglen != shared_siglen)break;
        if (sigarr[i]->signame != shared_signame)break;
    }
    if (i == sig->siglen) {
        output = (genosig_t*) malloc(sizeof(genosig_t));
        memcpy(output, sig, sizeof(genosig_t));
        output->sigdata = NULL;
        output->siglen = shared_siglen;
        output->sigtype = shared_sigtype;
        output->signame = shared_signame;
        newname[0] = "Average";
        newname[1] = sig->name;
        if (newname[1]) output->name = text_join(newname, ".", NULL, NULL, 2);
        else if (sig->filepath) {
            newname[1] = output->filepath;
            output->name = text_join(newname, ".", NULL, NULL, 2);
        }
        else output->name = text_join(newname, ".", NULL, NULL, 1);
        output->allocatedname = 1;
        tmpd = calloc(shared_siglen, sizeof(double));
        totweight = 0;
        switch (shared_sigtype) {
            case _SIGTYPE_DBL:
                for (i = 0;i < sig->siglen;i++) {
                    if (is_weighted) {
                        tmpd2 = vec_cpy((double*)(sigarr[i]->sigdata), shared_siglen);
                        vec_scale(tmpd2, sigarr[i]->weight, shared_siglen);
                        vec_add(tmpd, tmpd2, shared_siglen);
                        free(tmpd2);
                        totweight += sigarr[i]->weight;
                    }
                    else {
                        vec_add(tmpd, (double*)(sigarr[i]->sigdata), shared_siglen);
                        totweight += 1.0;
                    }
                }
                vec_scale(tmpd, 1 / totweight, shared_siglen);
                output->sigsize = shared_siglen * sizeof(double);
                output->sigdata = (void*)tmpd;
                break;
            case _SIGTYPE_I32:
                for (i = 0;i < sig->siglen;i++) {
                    tmpd2 = veci32_to_vec((int32_t*)(sigarr[i]->sigdata), shared_siglen);
                    if (is_weighted) {
                        vec_scale(tmpd2, sigarr[i]->weight, shared_siglen);
                        totweight += sigarr[i]->weight;
                    }
                    else {
                        totweight += 1.0;
                    }
                    vec_add(tmpd, tmpd2, shared_siglen);
                    free(tmpd2);
                }
                vec_scale(tmpd, 1 / totweight, shared_siglen);
                output->sigsize = shared_siglen * sizeof(int32_t);
                output->sigdata = (void*)vec_to_veci32(tmpd,shared_siglen);
                free(tmpd);
                break;
            case _SIGTYPE_I64:
                for (i = 0;i < sig->siglen;i++) {
                    tmpd2 = veci64_to_vec((int64_t*)(sigarr[i]->sigdata), shared_siglen);
                    if (is_weighted) {
                        vec_scale(tmpd2, sigarr[i]->weight, shared_siglen);
                        totweight += sigarr[i]->weight;
                    }
                    else {
                        totweight += 1.0;
                    }
                    vec_add(tmpd, tmpd2, shared_siglen);
                    free(tmpd2);
                }
                vec_scale(tmpd, 1 / totweight, shared_siglen);
                output->sigsize = shared_siglen * sizeof(int64_t);
                output->sigdata = (void*)vec_to_veci32(tmpd, shared_siglen);
                free(tmpd);
                break;
            default:
                free(tmpd);
                free(output);
                output = NULL;
                break;
        }
    }
    else
        output = NULL;
    return output;
}
genosig_t* genosig_varsig(genosig_t* sig, size_t is_weighted) {
    /* variance signature of a non-multiplexed signature is always NULL */
    /* a new signature is generated otehrwise*/
    int16_t shared_sigtype;
    int16_t shared_signame;
    size_t shared_siglen;
    size_t i,j;
    genosig_t** sigarr;
    genosig_t* output;
    double* tmpd;
    double* weightvec;
    double* varvec;
    double totweight;
    char* newname[2];
    if (sig->sigtype != _SIGTYPE_MUXED) return NULL;
    sigarr = sig->sigdata;
    shared_sigtype = sigarr[0]->sigtype;
    shared_siglen = sigarr[0]->siglen;
    shared_signame = sigarr[0]->signame;
    for (i = 0;i < sig->siglen;i++) {
        if (sigarr[i]->sigtype != shared_sigtype)break;
        if (sigarr[i]->siglen != shared_siglen)break;
        if (sigarr[i]->signame != shared_signame)break;
    }
    output = NULL;
    if (i == sig->siglen) {
        output = (genosig_t*)malloc(sizeof(genosig_t));
        memcpy(output, sig, sizeof(genosig_t));
        output->sigdata = NULL;
        output->siglen = shared_siglen;
        output->sigtype = shared_sigtype;
        output->signame = shared_signame;
        newname[0] = "EltVariance";
        newname[1] = sig->name;
        if (newname[1]) output->name = text_join(newname, ".", NULL, NULL, 2);
        else output->name = text_join(newname, ".", NULL, NULL, 1);
        output->allocatedname = 1;
        tmpd = (double*) malloc(shared_siglen*sig->siglen * sizeof(double));
        weightvec = (double*)malloc(sig->siglen * sizeof(double));
        varvec = (double*)malloc(shared_siglen * sizeof(double));
        totweight = 0;
        if (is_weighted) {
            for (j = 0;j < sig->siglen;j++) {
                weightvec[j] = sigarr[j]->weight;
            }
            totweight = vec_sum(weightvec, sig->siglen);
        }
        else {
            for (j = 0;j < sig->siglen;j++) {
                weightvec[j] = 1.0;
            }
            totweight = (double)(sig->siglen);
        }
        if (totweight == 0)totweight = 1;
        switch (shared_sigtype) {
        case _SIGTYPE_DBL:
            for (j = 0;j < sig->siglen;j++) {
                for (i = 0;i < shared_siglen;i++) {
                    tmpd[i*sig->siglen + j] = ((double*)(sigarr[j]->sigdata))[i];
                }
            }
            for (i = 0;i < shared_siglen;i++) {
                varvec[i] = vec_variance_weighted(tmpd + i*sig->siglen, weightvec, sig->siglen);
            }
            output->sigsize = shared_siglen * sizeof(double);
            output->sigdata = (void*)varvec;
            break;
        case _SIGTYPE_I32:
            for (j = 0;j < sig->siglen;j++) {
                for (i = 0;i < shared_siglen;i++) {
                    tmpd[i*sig->siglen + j] = (double)(((int32_t*)(sigarr[j]->sigdata))[i]);
                }
            }
            for (i = 0;i < shared_siglen;i++) {
                varvec[i] = vec_variance_weighted(tmpd + i*sig->siglen, weightvec, sig->siglen);
            }
            output->sigsize = shared_siglen * sizeof(int32_t);
            output->sigdata = (void*)vec_to_veci32(varvec, shared_siglen);
            free(varvec);
            break;
        case _SIGTYPE_I64:
            for (j = 0;j < sig->siglen;j++) {
                for (i = 0;i < shared_siglen;i++) {
                    tmpd[i*sig->siglen + j] = (double)(((int64_t*)(sigarr[j]->sigdata))[i]);
                }
            }
            for (i = 0;i < shared_siglen;i++) {
                varvec[i] = vec_variance_weighted(tmpd + i*sig->siglen, weightvec, sig->siglen);
            }
            output->sigsize = shared_siglen * sizeof(int64_t);
            output->sigdata = (void*)vec_to_veci64(varvec, shared_siglen);
            free(varvec);
            break;
        default:
            free(tmpd);
            free(output);
            output = NULL;
            break;
        }
        free(tmpd);
        free(weightvec);
    }
    else
        output = NULL;
    return output;
}
genosig_t* genosig_normparamsig(genosig_t* sig, size_t is_weighted) {
    /* calculates mu and sigma for each component of the signature */
    /* since the variance of a non-multiplexed signature is NULL, so is the output of this function */
    int16_t shared_sigtype;
    int16_t shared_signame;
    size_t shared_siglen;
    size_t i, j;
    genosig_t** sigarr;
    genosig_t* output;
    double* tmpd;
    double* tmpd2;
    double* varvec;
    double* weightvec;
    double totweight;
    char* newname[2];
    char* data_c;
    void* avgvec_v;
    void* varvec_v;
    size_t avgvec_size;
    size_t varvec_size;
    if (sig->sigtype != _SIGTYPE_MUXED) return NULL;
    sigarr = sig->sigdata;
    shared_sigtype = sigarr[0]->sigtype;
    shared_siglen = sigarr[0]->siglen;
    shared_signame = sigarr[0]->signame;
    for (i = 0;i < sig->siglen;i++) {
        if (sigarr[i]->sigtype != shared_sigtype)break;
        if (sigarr[i]->siglen != shared_siglen)break;
        if (sigarr[i]->signame != shared_signame)break;
    }
    if (i == sig->siglen) {
        output = (genosig_t*)malloc(sizeof(genosig_t));
        memcpy(output, sig, sizeof(genosig_t));
        output->sigdata = NULL;
        output->siglen = shared_siglen;
        output->sigtype = shared_sigtype;
        output->signame = _SIGNAME_MULTINORM;
        newname[0] = "NormParams";
        newname[1] = sig->name;
        if (newname[1]) output->name = text_join(newname, ".", NULL, NULL, 2);
        else if (sig->filepath) {
            newname[1] = output->filepath;
            output->name = text_join(newname, ".", NULL, NULL, 2);
        }
        else output->name = text_join(newname, ".", NULL, NULL, 1);
        output->allocatedname = 1;
        tmpd = calloc(shared_siglen, sizeof(double));
        weightvec = (double*)malloc(sig->siglen * sizeof(double));
        varvec = (double*)malloc(shared_siglen * sizeof(double));
        totweight = 0;
        if (is_weighted) {
            for (j = 0;j < sig->siglen;j++) {
                weightvec[j] = sigarr[j]->weight;
            }
            totweight = vec_sum(weightvec, sig->siglen);
        }
        else {
            for (j = 0;j < sig->siglen;j++) {
                weightvec[j] = 1.0;
            }
            totweight = (double)(sig->siglen);
        }
        switch (shared_sigtype) {
            case _SIGTYPE_DBL:
                for (i = 0;i < sig->siglen;i++) {
                    if (is_weighted) {
                        tmpd2 = vec_cpy((double*)(sigarr[i]->sigdata), shared_siglen);
                        vec_scale(tmpd2, sigarr[i]->weight, shared_siglen);
                        vec_add(tmpd, tmpd2, shared_siglen);
                        free(tmpd2);
                        totweight += sigarr[i]->weight;
                    }
                    else {
                        vec_add(tmpd, (double*)(sigarr[i]->sigdata), shared_siglen);
                        totweight += 1.0;
                    }
                }
                vec_scale(tmpd, 1 / totweight, shared_siglen);
                avgvec_size = shared_siglen * sizeof(double);
                avgvec_v = (void*)tmpd;
                break;
            case _SIGTYPE_I32:
                for (i = 0;i < sig->siglen;i++) {
                    tmpd2 = veci32_to_vec((int32_t*)(sigarr[i]->sigdata), shared_siglen);
                    if (is_weighted) {
                        vec_scale(tmpd2, sigarr[i]->weight, shared_siglen);
                        totweight += sigarr[i]->weight;
                    }
                    else {
                        totweight += 1.0;
                    }
                    vec_add(tmpd, tmpd2, shared_siglen);
                    free(tmpd2);
                }
                vec_scale(tmpd, 1 / totweight, shared_siglen);
                avgvec_size = shared_siglen * sizeof(int32_t);
                avgvec_v = (void*)vec_to_veci32(tmpd, shared_siglen);
                free(tmpd);
                break;
            case _SIGTYPE_I64:
                for (i = 0;i < sig->siglen;i++) {
                    tmpd2 = veci64_to_vec((int64_t*)(sigarr[i]->sigdata), shared_siglen);
                    if (is_weighted) {
                        vec_scale(tmpd2, sigarr[i]->weight, shared_siglen);
                        totweight += sigarr[i]->weight;
                    }
                    else {
                        totweight += 1.0;
                    }
                    vec_add(tmpd, tmpd2, shared_siglen);
                    free(tmpd2);
                }
                vec_scale(tmpd, 1 / totweight, shared_siglen);
                avgvec_size = shared_siglen * sizeof(int64_t);
                avgvec_v = (void*)vec_to_veci32(tmpd, shared_siglen);
                free(tmpd);
                break;
            default:
                avgvec_size = 0;
                avgvec_v = NULL;
                break;
        }
        tmpd = (double*)malloc(shared_siglen*sig->siglen * sizeof(double));
        switch (shared_sigtype) {
            case _SIGTYPE_DBL:
                for (j = 0;j < sig->siglen;j++) {
                    for (i = 0;i < shared_siglen;i++) {
                        tmpd[i*sig->siglen + j] = ((double*)(sigarr[j]->sigdata))[i];
                    }
                }
                for (i = 0;i < shared_siglen;i++) {
                    varvec[i] = sqrt(vec_variance_weighted(tmpd + i*sig->siglen, weightvec, sig->siglen));
                }
                free(tmpd);
                varvec_size = shared_siglen * sizeof(double);
                varvec_v = (void*)varvec;
                break;
            case _SIGTYPE_I32:
                for (j = 0;j < sig->siglen;j++) {
                    for (i = 0;i < shared_siglen;i++) {
                        tmpd[i*sig->siglen + j] = (double)(((int32_t*)(sigarr[j]->sigdata))[i]);
                    }
                }
                for (i = 0;i < shared_siglen;i++) {
                    varvec[i] = sqrt(vec_variance_weighted(tmpd + i*sig->siglen, weightvec, sig->siglen));
                }
                free(tmpd);
                varvec_size = shared_siglen * sizeof(int32_t);
                varvec_v = (void*)vec_to_veci32(varvec, shared_siglen);
                free(varvec);
                break;
            case _SIGTYPE_I64:
                for (j = 0;j < sig->siglen;j++) {
                    for (i = 0;i < shared_siglen;i++) {
                        tmpd[i*sig->siglen + j] = (double)(((int64_t*)(sigarr[j]->sigdata))[i]);
                    }
                }
                for (i = 0;i < shared_siglen;i++) {
                    varvec[i] = sqrt(vec_variance_weighted(tmpd + i*sig->siglen, weightvec, sig->siglen));
                }
                free(tmpd);
                varvec_size = shared_siglen * sizeof(int64_t);
                varvec_v = (void*)vec_to_veci64(varvec, shared_siglen);
                free(varvec);
                break;
            default:
                varvec_size = 0;
                varvec_v = NULL;
                break;
        }
        if (!avgvec_v || !varvec_v) {
            output->sigdata = NULL;
            output->sigsize = 0;
            genosig_free(output);
            output = NULL;
        }
        else {
            data_c = (char*)malloc(avgvec_size + varvec_size);
            memcpy(data_c, avgvec_v, avgvec_size);
            memcpy(data_c+avgvec_size, varvec_v, varvec_size);
            free(varvec_v);
            free(avgvec_v);
            output->sigsize = avgvec_size + varvec_size;
            output->sigdata = (void*)data_c;
        }
    }
    else
        output = NULL;
    return output;
}
static inline int64_t _genosig_nextrckmer(int64_t kmerrc, int64_t firstmask, int64_t firstneg) {
    int64_t mask, neg;
    mask = firstmask;
    neg = firstneg;
    while (mask && (kmerrc&mask) == 0) {
        kmerrc += mask;
        mask >>= 2;
        neg >>= 2;
    }
    kmerrc -= neg;
    return kmerrc;
}
genosig_t* genosig_keepnonredundant(genosig_t* sig, size_t k) {
    int64_t kmer,kmerrc,nkmers, j;
    int64_t firstmask,firstneg;
    int32_t* sigarri32;
    int64_t* sigarri64;
    double* sigarrdbl;
    nkmers = (1LL << (2 * k));
    if ((int64_t)(sig->siglen) == nkmers) {
        kmerrc = 0;
        firstneg = (nkmers >> 2);
        firstmask = nkmers - firstneg;
        kmerrc = nkmers - 1;
        j = 0;
        switch (sig->sigtype) {
        case _SIGTYPE_I32:
            sigarri32 = sig->sigdata;
            for (kmer = 0;kmer < nkmers;kmer++) {
                if (kmerrc >= kmer) {
                    sigarri32[j] = sigarri32[kmer];
                    j++;
                }
                kmerrc = _genosig_nextrckmer(kmerrc, firstmask, firstneg);
            }
            sig->siglen = j;
            sig->sigsize = j * sizeof(int32_t);
            break;
        case _SIGTYPE_I64:
            sigarri64 = sig->sigdata;
            for (kmer = 0;kmer < nkmers;kmer++) {
                if (kmerrc >= kmer) {
                    sigarri64[j] = sigarri64[kmer];
                    j++;
                }
                kmerrc = _genosig_nextrckmer(kmerrc, firstmask, firstneg);
            }
            sig->siglen = j;
            sig->sigsize = j * sizeof(int64_t);
            break;
        case _SIGTYPE_DBL:
            sigarrdbl = sig->sigdata;
            for (kmer = 0;kmer < nkmers;kmer++) {
                if (kmerrc >= kmer) {
                    sigarrdbl[j] = sigarrdbl[kmer];
                    j++;
                }
                kmerrc = _genosig_nextrckmer(kmerrc, firstmask, firstneg);
            }
            sig->siglen = j;
            sig->sigsize = j * sizeof(double);
            break;
        default:
            break;
        }
    }
    return sig;
}

/* distance */
static inline void _scale_vec_by_mode(double* vec, double* sigA, double* sigB, size_t len, int scalingmode) {
    double dlen;
    double tmp;
    dlen = (double)len;
    if (scalingmode == GENODIST_NORMELTSBYAVG) {
        tmp = 0.0;
        tmp += vec_avg(sigA, len);
        tmp += vec_avg(sigB, len);
        tmp /= 2.0;
        vec_scale(vec, 1.0 / tmp, len);
    }
    else if (scalingmode == GENODIST_NORMELTSBYDIMS) {
        vec_scale(vec, 1.0 / dlen, len);
    }
    else if (scalingmode == GENODIST_NORMELTSBYINVDIMS) {
        vec_scale(vec, dlen, len);
    }
}
static inline void _scale_veci32_by_mode(double* vec, int32_t* sigA, int32_t* sigB, size_t len, int scalingmode) {
    double dlen;
    double tmp;
    dlen = (double)len;
    if (scalingmode == GENODIST_NORMELTSBYAVG) {
        tmp = 0.0;
        tmp += (double)veci32_avg(sigA, len);
        tmp += (double)veci32_avg(sigB, len);
        tmp /= 2.0;
        vec_scale(vec, 1.0 / tmp, len);
    }
    else if (scalingmode == GENODIST_NORMELTSBYDIMS) {
        vec_scale(vec, 1.0 / dlen, len);
    }
    else if (scalingmode == GENODIST_NORMELTSBYINVDIMS) {
        vec_scale(vec, dlen, len);
    }
}
double genodist_bytaxonomy(genosig_t* A, genosig_t* B, any_t flags) {
    if (A->lineage && B->lineage) return lineage_majornamerelatedness(A->lineage, B->lineage);
    return -1.0;
}
double genodist_manhattan(genosig_t* A, genosig_t* B, any_t flags){
    double result;
    void* diff;
    double* diff_d;
    
    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;

    if (A->sigtype == _SIGTYPE_I32) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;
        diff = malloc(A->sigsize);
        memcpy(diff, A->sigdata, A->sigsize);

        veci32_subtract((int32_t*)diff, (int32_t*)(B->sigdata), A->siglen);
        diff_d = veci32_to_vec(diff, A->siglen);
        free(diff);
        _scale_veci32_by_mode(diff_d, (int32_t*)(A->sigdata), (int32_t*)(B->sigdata), A->siglen, flags.u8[0]);
        vec_abs(diff_d, A->siglen);
        result = vec_sum(diff_d, A->siglen);
        free(diff_d);
    }
    else if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;

        diff_d = vec_cpy(A->sigdata, A->siglen);
        vec_subtract(diff_d, (double*)(B->sigdata), A->siglen);

        _scale_vec_by_mode(diff_d, (double*)(A->sigdata), (double*)(B->sigdata), A->siglen, flags.u8[0]);
        vec_abs(diff_d, A->siglen);
        result = vec_sum(diff_d, A->siglen);
        free(diff_d);
    }
    return result;
}
double genodist_euclidian(genosig_t* A, genosig_t* B, any_t flags){
    double result;
    void* diff;
    double* diff_d;

    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;

    if (A->sigtype == _SIGTYPE_I32) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;

        diff = malloc(A->sigsize);
        memcpy(diff, A->sigdata, A->sigsize);
        veci32_subtract((int32_t*)diff, (int32_t*)(B->sigdata), A->siglen);
        diff_d = veci32_to_vec(diff, A->siglen);
        free(diff);

        _scale_veci32_by_mode(diff_d, (int32_t*)(A->sigdata), (int32_t*)(B->sigdata), A->siglen, flags.u8[0]);
        vec_dot(diff_d, diff_d, A->siglen);
        result = vec_sum(diff_d, A->siglen);
        free(diff_d);
    }
    else if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;

        diff_d = vec_cpy(A->sigdata, A->siglen);
        vec_subtract(diff_d, (double*)(B->sigdata), A->siglen);
        
        _scale_vec_by_mode(diff_d, (double*)(A->sigdata), (double*)(B->sigdata), A->siglen, flags.u8[0]);
        vec_dot(diff_d, diff_d, A->siglen);
        result = vec_sum(diff_d, A->siglen);
        free(diff_d);
    }
    return result;
}
double genodist_pearscorr(genosig_t* A, genosig_t* B, any_t flags){
    double result;
    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;

    if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;
        result = vec_pearsoncorr((double*)(A->sigdata), (double*)(B->sigdata), A->siglen);
    }
    /* this ensures that the output range is [0:1] with identical genomes having a "distance" of 0 */
    return 1.0 - (result + 1.0) / 2.0;
}
double genodist_pearscorr_unbound(genosig_t* A, genosig_t* B, any_t flags) {
    double result;
    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;

    if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;
        result = vec_pearsoncorr((double*)(A->sigdata), (double*)(B->sigdata), A->siglen);
    }
     return result;
}
double genodist_rankcorr(genosig_t* A, genosig_t* B, any_t flags){
    double result;
    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;

    if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;
        result = vec_spearmancorr((double*)(A->sigdata), (double*)(B->sigdata), A->siglen);
    }
    /* this ensures that the output range is [0:1] with identical genomes having a "distance" of 0 */
    return 1.0 - (result + 1.0) / 2.0;
}
double genodist_rankcorr_unbound(genosig_t* A, genosig_t* B, any_t flags) {
    double result;
    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;

    if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;
        result = vec_spearmancorr((double*)(A->sigdata), (double*)(B->sigdata), A->siglen);
    }
    return result;
}
double genodist_jaccard(genosig_t* A, genosig_t* B, any_t flags){
    return -1;
}
double genodist_hamming(genosig_t* A, genosig_t* B, any_t flags) {
    size_t diffs;
    size_t i;
    size_t eltsize;
    char* zeroelt;
    char* sigA;
    char* sigB;
    size_t nonzeroA;
    size_t nonzeroB;
    double result;

    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return -1;

    eltsize = A->sigsize / A->siglen;

    diffs = 0;
    sigA = (char*)A->sigdata;
    sigB = (char*)B->sigdata;
    for (i = 0;i < A->siglen;i++) {
        if (memcmp(sigA+i*eltsize, sigB + i*eltsize, eltsize))diffs++;
    }
    result = (double)diffs;
    if (flags.u8[0] == GENODIST_NORMELTSBYDIMS) result /= (double)(A->siglen);
    else if (flags.u8[0] == GENODIST_NORMELTSBYINVDIMS) result *= (double)(A->siglen);
    else if (flags.u8[0] == GENODIST_NORMELTSBYAVG) {
        zeroelt = calloc(eltsize,1);
        nonzeroA = 0;
        nonzeroB = 0;
        if (A->sigtype == _SIGTYPE_DBL) {
            /* IEEE 754 allows two representations of 0, hence the separation */
            for (i = 0;i < A->siglen;i++) {
                if (*((double*)(sigA + i*eltsize)) != 0.0) nonzeroA++;
            }
            for (i = 0;i < B->siglen;i++) {
                if (*((double*)(sigB + i*eltsize)) != 0.0) nonzeroB++;
            }
        }
        else {
            /* for all non-floating point formats, we assume that only 0 <=> nothing */
            for (i = 0;i < A->siglen;i++) {
                if (memcmp(sigA + i*eltsize, zeroelt, eltsize))nonzeroA++;
            }
            for (i = 0;i < B->siglen;i++) {
                if (memcmp(sigB + i*eltsize, zeroelt, eltsize))nonzeroB++;
            }
        }
        free(zeroelt);
        result = 2.0 * result / ((double)(nonzeroA + nonzeroB));
    }
    return result;
}
static inline double resemblance2ANI(double v, double k, double genomesize) {
    if (v == 0.0)return 0;
    return 1.0 + log(v * 2 / (1 + v)) / k;
}
double genodist_approxANI(genosig_t* A, genosig_t* B, any_t flags){
    int64_t* s1;
    int64_t* s2;
    size_t i, j;
    size_t AinterB;
    size_t AunionB;
    double k;
    double resemblance;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return -1;

    k = flags.d;
    if (k == 0) {
        if (A->window <= MINMINHASHLEN)k = 21;
        else k = (double)(A->window);
    }
    s1 = (int64_t*)A->sigdata;
    s2 = (int64_t*)B->sigdata;
    i = j = 0;
    AinterB = AunionB = 0;
    while (i < A->siglen && j < B->siglen) {
        if (s1[i] == s2[j]) {
            i++;
            j++;
            AinterB++;
        }
        else if (s1[i] > s2[j]) {
            j++;
        }
        else {
            i++;
        }
        AunionB++;
    }
    resemblance = ((double)AinterB) / ((double)AunionB);
    return 1 - resemblance2ANI(resemblance, k, (double)((s1[0] + s2[0]) / 2));
}
double genodist_approxANI_unbound(genosig_t* A, genosig_t* B, any_t flags) {
    double tmp;
    tmp = 100. - (genodist_approxANI(A, B, flags) * 100);
    if (tmp > 100)tmp = 0;
    return tmp;
}
double genodist_SVC(genosig_t* A, genosig_t* B, any_t maxdelta) {
    /* similar value count is a count of values where diff < maxdelta */
    double result;
    void* diff;
    double* diff_d;
    size_t i;
    size_t count;
    size_t siglen_to_consider;

    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;
    siglen_to_consider = A->siglen;
    if (A->signame == _SIGNAME_KARLINL) siglen_to_consider--;

    if (A->sigtype == _SIGTYPE_I32) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;
        diff = malloc(A->sigsize);
        memcpy(diff, A->sigdata, A->sigsize);

        veci32_subtract((int32_t*)diff, (int32_t*)(B->sigdata), siglen_to_consider);
        diff_d = veci32_to_vec(diff, siglen_to_consider);
        free(diff);
        vec_abs(diff_d, siglen_to_consider);
        count = 0;
        for (i = 0; i < siglen_to_consider;i++) {
            if (diff_d[i]>(double)(maxdelta.i64)) count++;
        }
        result = (double)count/(double)(siglen_to_consider);
        free(diff_d);
    }
    else if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;

        diff_d = vec_cpy(A->sigdata, siglen_to_consider);
        vec_subtract(diff_d, (double*)(B->sigdata), siglen_to_consider);

        vec_abs(diff_d, siglen_to_consider);
        count = 0;
        for (i = 0; i < siglen_to_consider;i++) {
            if (diff_d[i]>(double)(maxdelta.d)) count++;
        }
        result = (double)count / (double)(siglen_to_consider);
        free(diff_d);
    }
    return result;
}
double genodist_satman(genosig_t* A, genosig_t* B, any_t maxdelta){
    /*
    This "saturated manhattan distance" was derived as a countinuous
    variant of similar value count above, and seemed to yield better results
    */
    double result;
    void* diff;
    double* diff_d;
    size_t i;
    size_t maxed;
    size_t siglen_to_consider;

    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;
    siglen_to_consider = A->siglen;
    if (A->signame == _SIGNAME_KARLINL) siglen_to_consider--;

    maxed = 0;
    if (A->sigtype == _SIGTYPE_I32) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;
        diff = malloc(A->sigsize);
        memcpy(diff, A->sigdata, A->sigsize);

        veci32_subtract((int32_t*)diff, (int32_t*)(B->sigdata), siglen_to_consider);
        diff_d = veci32_to_vec(diff, siglen_to_consider);
        free(diff);
        vec_abs(diff_d, siglen_to_consider);
        for (i = 0; i < siglen_to_consider;i++) {
            if (diff_d[i] > (double)(maxdelta.i64)) {
                diff_d[i] = (double)(maxdelta.i64);
                maxed++;
            }
        }
        result = vec_sum(diff_d, siglen_to_consider);
        free(diff_d);
    }
    else if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;

        diff_d = vec_cpy(A->sigdata, siglen_to_consider);
        vec_subtract(diff_d, (double*)(B->sigdata), siglen_to_consider);

        vec_abs(diff_d, siglen_to_consider);
        for (i = 0; i < siglen_to_consider;i++) {
            if (diff_d[i] > maxdelta.d) {
                diff_d[i] = maxdelta.d;
                maxed++;
            }
        }
        result = vec_sum(diff_d, siglen_to_consider);
        free(diff_d);
    }
    if (maxdelta.d > 0) {
        result /= (double)(siglen_to_consider);
        result /= maxdelta.d;
    }
    else if (maxdelta.d == 0) {
        result = (double)maxed;
        result /= (double)(siglen_to_consider);
    }
    return result;
}
double genodist_PaSiTL(genosig_t* A, genosig_t* B, any_t maxdelta) {
    /*
    This "saturated manhattan distance" was derived as a countinuous
    variant of similar value count above, and seemed to yield better results
    */
    double result;
    double* diff_d;
    double* Asig;
    double* Bsig;
    double Alen;
    double Blen;
    double avglen;
    double maxdiff;
    size_t i;
    size_t maxed;
    size_t siglen_to_consider;

    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype || A->sigtype != _SIGTYPE_DBL || A->signame != _SIGNAME_KARLINL) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;
    siglen_to_consider = A->siglen-1;
    
    maxed = 0;
    /* only compare vectors with the same length */
    if (A->siglen != B->siglen) return result;

    Asig = vec_cpy((double*)(A->sigdata), A->siglen);
    Bsig = vec_cpy((double*)(B->sigdata), B->siglen);
    Alen = Asig[siglen_to_consider];
    Blen = Bsig[siglen_to_consider];
    avglen = (Alen + Blen) / 2;
    
    for (i = 0; i < siglen_to_consider;i++) {
        Asig[i] = log(Asig[i]);
        Bsig[i] = log(Bsig[i]);
    }
    vec_scale(Asig, log(Alen), siglen_to_consider);
    vec_scale(Bsig, log(Blen), siglen_to_consider);
    
    diff_d = Asig;
    vec_subtract(diff_d, Bsig, siglen_to_consider);
    free(Bsig);
    vec_abs(diff_d, siglen_to_consider);

    maxdiff = log(1+maxdelta.d)*log(avglen);
    for (i = 0; i < siglen_to_consider;i++) {
        if (diff_d[i] > maxdiff) {
            diff_d[i] = maxdiff;
            maxed++;
        }
    }
    result = vec_sum(diff_d, siglen_to_consider);
    free(diff_d);
    if (maxdelta.d > 0) {
        result /= (double)(siglen_to_consider);
        result /= maxdiff;
    }
    else if (maxdelta.d == 0) {
        result = (double)maxed;
        result /= (double)(siglen_to_consider);
    }

    return result;
}
double genodist_sateucl(genosig_t* A, genosig_t* B, any_t maxdelta) {
    /*
    This "saturated manhattan distance" was derived as a countinuous
    variant of similar value count above, and seemed to yield better results
    */
    double result;
    void* diff;
    double* diff_d;
    double mx2;
    size_t i;
    size_t maxed;
    size_t siglen_to_consider;

    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (A->sigtype != B->sigtype) return result;
    if (A->siglen == 0 || B->siglen == 0) return result;

    siglen_to_consider = A->siglen;
    if (A->signame == _SIGNAME_KARLINL) siglen_to_consider--;

    maxed = 0;
    mx2 = maxdelta.d*maxdelta.d;
    if (A->sigtype == _SIGTYPE_I32) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;
        diff = malloc(A->sigsize);
        memcpy(diff, A->sigdata, A->sigsize);

        veci32_subtract((int32_t*)diff, (int32_t*)(B->sigdata), siglen_to_consider);
        diff_d = veci32_to_vec(diff, siglen_to_consider);
        free(diff);
        vec_abs(diff_d, siglen_to_consider);
        for (i = 0; i < siglen_to_consider;i++) {
            if (diff_d[i] > (double)(maxdelta.i64)) {
                diff_d[i] = (double)(maxdelta.i64);
                maxed++;
            }
        }
        result = vec_sum(diff_d, siglen_to_consider);
        free(diff_d);
    }
    else if (A->sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (A->siglen != B->siglen) return result;

        diff_d = vec_cpy(A->sigdata, siglen_to_consider);
        vec_subtract(diff_d, (double*)(B->sigdata), siglen_to_consider);

        vec_dot(diff_d, diff_d, siglen_to_consider);
        for (i = 0; i < siglen_to_consider;i++) {
            if (diff_d[i] > mx2) {
                diff_d[i] = mx2;
                maxed++;
            }
        }
        result = vec_sum(diff_d, siglen_to_consider);
        free(diff_d);
    }
    if (maxdelta.d > 0) {
        result /= (double)(siglen_to_consider);
        result /= mx2;
    }
    else if (maxdelta.d == 0) {
        result = (double)maxed;
        result /= (double)(siglen_to_consider);
    }
    return sqrt(result);
}
double _genodist_NI(suftree_t* s1, suftree_t* s2, size_t* alnlen) {
    char* seq2;
    size_t seq2len;
    char* seq1;
    size_t seq1len;
    if (!s1 || !s2) {
        *alnlen = 1;
        return 1;
    }
    seq2 = suftree_getstrptr(s2, &seq2len);
    seq1 = suftree_getstrptr(s1, &seq1len);
    if (seq1len < 100 || seq1len>seq2len * 2) return 1.0;
    if (seq2len < 100 || seq2len>seq1len * 2) return 1.0;
    return suftree_roughalign(seq2, seq2len, s1, 4+(size_t)log10((double)(seq1len+seq2len)), 1, 1, alnlen);
}
double genodist_ANI(genosig_t* A, genosig_t* B, any_t flags){
    size_t i;
    suftree_t* stA;
    suftree_t* stB;
    double ANI;
    double NI;
    size_t totlen;
    size_t alnlen;
    if (A->sigtype != B->sigtype || A->sigtype != _SIGTYPE_SUFTREEARRAY) return -1.0;
    totlen = 0;
    ANI = 0;
    for (i = 0;i < A->siglen;i++) {
        stA = ((suftree_t**)(A->sigdata))[i];
        stB = ((suftree_t**)(B->sigdata))[i];
        NI = _genodist_NI(stA, stB, &alnlen);
        if (flags.i8[0] == GENODIST_NORMBYSEQ) ANI += NI;
        else ANI += NI * (double)alnlen;
        totlen += alnlen;
    }
    if (flags.i8[0] == GENODIST_NORMBYSEQ) {
        ANI /= (double)A->siglen;
    }
    else {
        ANI /= (double)totlen;
    }
    return ANI;
}
double genodist_samespecies(genosig_t* A, genosig_t* B, any_t flags){
    return -1;
}
double genodist_sametype(genosig_t* A, genosig_t* B, any_t flags){
    return -1;
}
double genodist_avgzscore(genosig_t* sample, genosig_t* expected_distribution, any_t unused) {
    double result;
    void* diff;
    double* diff_d;
    double* scaling_factor;
    genosig_t* tmpptr;
    int sigtype;
    size_t siglen;

    result = -1.0;
    /* if different signature types are supplied, don't bother comparing them */
    if (sample->sigtype != expected_distribution->sigtype) return result;

    /* if "expected_distribution" does actually represent distribution parameters, also don't bother doing anything else */
    if (sample->signame == _SIGNAME_MULTINORM && expected_distribution->signame != _SIGNAME_MULTINORM) {
        tmpptr = sample;
        sample = expected_distribution;
        expected_distribution = tmpptr;
    }
    if (expected_distribution->signame != _SIGNAME_MULTINORM) return result;

    sigtype = sample->sigtype;
    siglen = sample->siglen;
    if (sample->sigtype == _SIGTYPE_I32) {
        /* only compare vectors with the same length */
        if (siglen != expected_distribution->siglen) return result;
        diff = malloc(sample->sigsize);
        memcpy(diff, sample->sigdata, sample->sigsize);

        veci32_subtract((int32_t*)diff, (int32_t*)(expected_distribution->sigdata),siglen);
        diff_d = veci32_to_vec(diff, siglen);
        free(diff);
        vec_abs(diff_d, siglen);
        scaling_factor = veci32_to_vec(((int32_t*)(expected_distribution->sigdata)) + siglen, siglen);
        vec_inverse(scaling_factor, siglen);
        vec_dot(diff_d, scaling_factor, siglen);
        result = vec_avg(diff_d, siglen);
        free(diff_d);
        free(scaling_factor);
    }
    else if (sigtype == _SIGTYPE_DBL) {
        /* only compare vectors with the same length */
        if (siglen != expected_distribution->siglen) return result;

        diff_d = vec_cpy(sample->sigdata, siglen);
        vec_subtract(diff_d, (double*)(expected_distribution->sigdata), siglen);
        vec_abs(diff_d, siglen);
        scaling_factor = vec_cpy(((double*)(expected_distribution->sigdata)) + siglen, siglen);
        vec_inverse(scaling_factor, siglen);
        vec_dot(diff_d, scaling_factor, siglen);
        result = vec_avg(diff_d, siglen);
        free(diff_d);
        free(scaling_factor);
    }
    return result;
}

double genodist_externANIb_oneway(genosig_t* A, genosig_t* B, any_t blastpath){
    char* args[15] = { NULL,"-query",NULL,"-subject",NULL,"-outfmt","6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore","-perc_identity","50.0","-max_target_seqs","1000","-num_threads","1","-culling_limit","1" };
    char* resultstr;
    double resultdist;
    double NI;
    uint64_t totlen;
    uint64_t alnlen;
    size_t i,j;
    size_t reslen;
    int tries_remaining;
    args[0] = blastpath.str;
    args[2] = A->filepath;
    args[4] = B->filepath;
    tries_remaining = 3;
    resultstr = NULL;
    while(tries_remaining > 0 && resultstr == NULL) {
      resultstr = os_stdoutfromexec(blastpath.str, 15, args, &reslen, 0x10000);
      tries_remaining--;
    }
    if(!resultstr) return -1;
    i = j = 0;
    totlen = 0;
    resultdist = 0.0;
    while (i < reslen) {
        /* skip the first two fields */
        while (i < reslen-1 && resultstr[i] != '\t')i++;
        i++;
        if (i >= reslen - 1)break;
        while (i < reslen-1 && resultstr[i] != '\t')i++;
        i++;
        if (i >= reslen - 1)break;
        /* isolate the third field (perc_identity) */
        j = i;
        while (j < reslen-1 && resultstr[j] != '\t')j++;
        if (j >= reslen - 1)break;
        resultstr[j] = 0;
        NI = atof(resultstr + i);
        i = j+1;
        /* isolate the fourth field (aln_len) */
        while (j < reslen - 1 && resultstr[j] != '\t')j++;
        if (j >= reslen - 1)break;
        resultstr[j] = 0;
        alnlen = atoll(resultstr + i);
        while (i < reslen && resultstr[i] != '\n')i++;
        resultdist += (NI / 100.0)*((double)alnlen);
        totlen += alnlen;
    }
    free(resultstr);
    if(totlen==0) return -1;
    return ((double)totlen - resultdist) / (double)totlen;
}

double genodist_externANIb(genosig_t* A, genosig_t* B, any_t blastpath){
    return (genodist_externANIb_oneway(A,B,blastpath)+genodist_externANIb_oneway(B,A,blastpath))/2;
}
double genodist_externANIb_unbound(genosig_t* A, genosig_t* B, any_t blastpath) {
    double tmp;
    tmp = 100.-genodist_externANIb(A, B, blastpath)*100;
    if (tmp > 100)tmp = 0;
    return tmp;
}


double genodist_externdist(genosig_t* A, genosig_t* B, any_t progstring){
    static int errordone = 0;
    if (!errordone) {
        fprintf(stderr, "Please implement genodist_externdist for this to work!\n");
        errordone = 1;
    }
    return -1;
}
