// GLOBAL DEFINES
#define VERSION_NAME    "GenDisCal v1.0"
#ifdef _DEBUG
/*#define NOOPENMP*/
#endif
// INCLUDES
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#ifdef NOOPENMP
#define omp_get_max_threads()   1
#define omp_get_thread_num() 0
#define omp_set_num_threads(numthreads)
#else
#include <omp.h>
#endif

#include "datamap.h"
#include "osportstd.h"
#include "nucseq.h"
#include "vecops.h"
#include "argparser.h"
#include "GenDisCal_distances.h"
#include "textparsing.h"
#include "taxextractor.h"
#include "datatable.h"

int autothreads() {
    int mthreads;
    int mpow;
    mthreads = omp_get_max_threads();
    args_report_info(NULL, " (openmp) Max threads: %d\n", mthreads);
    if (mthreads <= 1) return 1;
    if (mthreads < 6)
        return mthreads - 1;
    if (mthreads < 16)
        return 4;
    mpow = 0;
    while (mthreads > 1) {
        mthreads >>= 1;
        mpow++;
    }
    return 1 << (mpow - 1);
}

args_t* GenDisCal_init_args(int argc, char** argv) {
    args_t* result;
    
    result = args_alloc();

    args_add(result, "basis", 'b', "str,int");
    args_add(result, "method", 'm', "str,float");
    args_add(result, "preset", 'p', "str");
    args_add(result, "filelist", 'f', "str");
    args_add(result, "search", 's', "str");
    args_add(result, "output", 'o', "str");
    args_add(result, "taxonomy", 't', "str,str");
    args_add(result, "histogram", 'g', "float");
    args_add(result, "distancematrix", 'x', "");
    args_add(result, "under", 'u', "float");
    args_add(result, "above", 'v', "float");
    args_add(result, "only", 'y', "str");

#ifndef NOOPENMP
    args_add(result, "nprocs", 'n', "int");
#endif

    args_add_help(result, "basis", "SIGNATURE BASIS",
        "Signature bases are used to transform a genome into a vector of numbers. "
        "Available signatures are as follows:\n"
        "   GC         -- use GC content as a signature\n"
        "   freq <k>   -- use frequencies of <k>-mers as a signature\n"
        "   karl <k>   -- use karlin signatures of <k>-mers\n"
        "   exp <k>    -- use <k>-mer bias with regards to monomer composition\n"
        /*"   bias <k>   -- use <k>-mer bias with regards to <k-1>-mer composition\n"*/
        "For single-stranded DNA and RNA, alternative versions of these functions are provided, "
        "beginning with 's_', as follows:\n"
        "   s_freq, s_karl, s_exp, s_bias\n",
        "signature basis to be used");
    args_add_help(result, "method", "DISTANCE COMPUTATION METHOD",
        "Distance computation methods are used to compute the dissimilarity between genome signatures."
        "Available methods are as follows:\n"
        "   manhattan/AMD  -- use the manhattan distance\n"
        "   euclid/ED      -- use the euclidian distance\n"
        "   hamming        -- use the hamming distance\n"
        "   PaSiT/SVC <f>  -- use the PaSiT distance with a non-standard threshold\n"
        "   pearson/corr   -- use Pearson correlation as dissimilarity measure\n",
        "distance calculation method to be used");
    args_add_help(result, "preset", "OPTION PRESET",
        "For convenience some option presets are provided, which have shown good preformance during "
        "the development and testing of this software, as well more complex methods which cannot easily be "
        "input normally. The following presets are provided:\n\n"
        " > PaSiT4 (-b karl 4 -m PaSiT 0.02)\n"
        "    The PaSiT4 method is suitable to conclude whether two genomes belong\n"
        "    to the same species.\n"
        "    The ambiguity region for the results is [0.22-0.36]\n"
        " > PaSiT6 (-b karl 6 -m PaSiT 0.02)\n"
        "    The PaSiT6 method is suitable to conclude whether two genomes belong\n"
        "    to the same type.\n"
        "    The ambiguity region for the results is [0.22-0.28]\n"
        " > TETRA (-b bias 4 -m corr)\n"
        "    TETRA is suitable to conclude whether two genomes belong to the same\n"
        "    species.\n"
        "    The ambiguity region for the results is [0.0005-0.005]\n"
 /*       " > program:<file> (no equivalent)\n"
        "    This will run an external program and format the result based on the remaining arguments."
        "    <file> should be formatted as follows:\n"
        "      line 1 | command line to summarize sequences (leave blank to use raw sequences)\n"
        "      line 2 | command line to compare summaries (or raw sequences)\n"
        "      line 3 | command line to transform the results into a tab-separated table\n"
        "    Some examples are provided, but need to be adapted manually to work:\n"
        "      <install path>/16S.eem\n"
        "      <install path>/MinHash.eem\n"
        "      <install path>/ANIb.eem\n"
        "    The software required to run some of these examples needs to be installed separately.\n"*/,
        "use a special preset of options (replaces -b and -m)");
    args_add_help(result, "filelist", "FILE LIST",
        VERSION_NAME " accepts genome files as free arguments, but when a large number of genomes "
        "need to be compared, it is often more useful to store a list of these genomes in a separate "
        "file. For " VERSION_NAME ", this file should contain one relative or absolute path to a "
        "readable (i.e. signature or nucleotide fasta) file per line.\n",
        "file containing a list of paths to files to use");
    args_add_help(result, "search", "SEARCH",
        "If present, this option will cause the program to compare the specified file against all "
        "other files, instead of performing all-against-all comparisons. If an invalid or no "
        "arguemnt is specified, the last indicated file will be considered to be the query. "
        "This is option disables replaces --distance matrix with a regular output.\n",
        "perform 1/all comparisons instead of all/all");
    args_add_help(result, "output", "OUTPUT FILE NAME",
        "Specifies the output file name.",
        "specify the output file name");
    args_add_help(result, "taxonomy", "TAXONOMY FILE",
        "Specifies the taxonomy file. This option is designed for testing whether the distance "
        "matches the expected taxonomical relation. The first argument should be the path to the "
        "taxonomy file, while the second should be the taxonomy format. Possible formats are:\n"
        "  PCOF_S  -- filename,PhylumClassOrderFamily_Genus_species_etc (default)\n"
        "  semicol -- filename;Phylum;Class;Order;Family;Genus;species;etc\n"
        "  pipe    -- filename | Phylum | Class | Order | Family | Genus | Species | etc\n",
        "taxonomy file");
    args_add_help(result, "histogram", "HISTOGRAM OUTPUT",
        "If this option is present, the output will be a multi-column histogram with the specified "
        "bin width (default: 0.001) instead of a list of comparisons. Rows represent bins, whereas "
        "columns represent the taxonomic relation based on the --taxonomy option.\n"
        "The output format is a always a comma-separated file.\n",
        "enable histogram output (overrides --distancematrix)");
    args_add_help(result, "distancematrix", "DISTANCE MATRIX OUTPUT",
        "If this option is present, the output will be a distance matrix instead of a list of "
        "comparisons.\n"
        "The output format is a always a comma-separated file.\n",
        "enable distance matrix output");
    args_add_help(result, "under", "VALUES UNDER A GIVEN VALUE",
        "Only report values below the specified value. (default:infinity)\n",
        "Only report values below the specified value");
    args_add_help(result, "above", "VALUES ABOVE A GIVEN VALUE",
        "Only report values above the specified value. (default:infinity)\n",
        "Only report values above the specified value");
    args_add_help(result, "only", "VALUES WITH A GIVEN TAXONOMIC RELATION",
        "If this option is present, only values associated with a given taxonomic relation "
        "will be reported. The available relations are as follows, and should be "
        "self-explanatory:\n"
        "  ?\n"
        "  Different_Phyla\n"
        "  Same_Phylum\n"
        "  Same_Class\n"
        "  Same_Order\n"
        "  Same_Family\n"
        "  Same_Genus\n"
        "  Same_Species\n"
        "  Same_Subspecies\n"
        "  Replicate\n",
        "Only report values with a given taxonomic relation");
#ifndef NOOPENMP
    args_add_help(result, "nprocs", "NUMBER OF PROCESSOR CORES TO USE",
        "This option defines the number of cores to use. By default, the number of cores "
        "used depends on how many are available:\n"
        "  1,2  :  1"
        "  3,4  :  All but one\n"
        "  5-15 :  4\n"
        "  16+  :  Half of the available cores, rounded down to a power of 2\n",
        "number of cores to use");
#endif
    args_parse(result, argc, argv);
    return result;
}
char** GenDisCal_get_file_list(args_t* args, size_t* nfiles) {
    size_t nf;
    size_t nl;
    size_t i, maxi;
    char* line;
    char** result;
    PF_t* f;
    nf = 0;
    result = NULL;
    if (args_ispresent(args, "filelist")) {
        line = args_getstr(args, "filelist", 0, NULL);
        if (!line) {
            args_report_info(NULL, "Argument for --filelist not specied, defaulting to \"files.list\"\n", line);
            line = "files.list";
        }
        PFopen(&f, line, "r");
        if (!f)args_report_warning(NULL, "File list <%s> could not be opened\n",line);
        else {
            nl = PFnumlines(f, 0);
            result = (char**)realloc(result, sizeof(char*)*nl);
            while (line = PFreadline(f)) {
                if (os_fileexists(line)) {
                    result[nf] = line;
                    nf++;
                }
                else {
                    if(strlen(line)>0) args_report_warning(NULL, "File <%s> does not exist\n", line);
                    free(line);
                }
            }
        }
    }
    nl = nf;
    maxi = args_countfreeargs(args);
    nl += maxi;
    result = (char**)realloc(result, sizeof(char*)*nl);
    for (i = 0;i < maxi;i++) {
        line = args_getstr(args, NULL, (int)i, NULL);
        if(line && os_fileexists(line)) {
            result[nf] = malloc(strlen(line)+1);
            memcpy(result[nf], line, strlen(line) + 1);
            nf++;
        }
        else {
            if (strlen(line)>0) args_report_warning(NULL, "File <%s> does not exist\n", line);
        }
    }
    if (args_ispresent(args, "search")) {
        line = args_getstr(args, "search", 0, NULL);
        if (line && os_fileexists(line)) {
            nl = nf + 1;
            result = (char**)realloc(result, sizeof(char*)*nl);
            result[nf] = malloc(strlen(line) + 1);
            memcpy(result[nf], line, strlen(line) + 1);
            nf++;
        }
    }
    *nfiles = nf;
    if(nf!=1)
        args_report_info(args,  _LLD_ " files were loaded.\n", nf);
    else
        args_report_info(args, "1 file was loaded.\n");
    return result;
}

typedef struct signature_t {
    void* data;
    size_t len;
    char* fname;
    int64_t lineage[9];
    int freefname;
} signature_t;
#define LINEAGE_PHYLUM      0
#define LINEAGE_CLASS       1
#define LINEAGE_ORDER       2
#define LINEAGE_FAMILY      3
#define LINEAGE_GENUS       4
#define LINEAGE_SPECIES     5
#define LINEAGE_SUBSPECIES  6
#define LINEAGE_STRAIN      7
#define LINEAGE_UNIQUEID    8

signature_t* signature_alloc() {
    signature_t* result;
    result = calloc(1, sizeof(signature_t));
    return result;
}
void signature_free(signature_t* target) {
    if (target) {
        if (target->data)free(target->data);
        if (target->fname && target->freefname) free(target->fname);
        target->data = NULL;
        target->fname = NULL;
        target->len = 0;
        free(target);
    }
}
int signature_taxsimlevel(signature_t* A, signature_t* B) {
    int i;
    i = 0;
    if (!A || !B || A->lineage[i] == 0 || B->lineage[i] == 0) return -1;
    while (i < 9 && A->lineage[i] == B->lineage[i] && A->lineage[i] != 0)
        i++;
    return i;
}

int set_basis_function(args_t* args, basis_function* bf, size_t* sig_len, size_t* kmer_len) {
    size_t maxi, i;
    size_t kmerlen;
    int nocalc;
    char* basisstr;
    char* presetstr;
    char* possible_bases[] = { "GC","freq","karl","exp","s_freq","s_karl","s_exp" };
    basis_function bf_list[] = { freqs,freqs,karsn,markz2,freqn,karln,markz1 };
    nocalc = 0;
    presetstr = args_getstr(args, "preset", 0, NULL);
    if (presetstr) {
        if (strcmp(presetstr, "TETRA") == 0) {
            *bf = TETRA;
            kmerlen = 4;
            *sig_len = (1LL << (4 * 2));
        }
        else if (strcmp(presetstr, "PaSiT4") == 0) {
            *bf = karsn;
            kmerlen = 4;
            *sig_len = (1LL << (4 * 2));
        }
        else if (strcmp(presetstr, "PaSiT6") == 0) {
            *bf = karsn;
            kmerlen = 6;
            *sig_len = (1LL << (6 * 2));
        }
        else {
            args_report_warning(NULL, "Unknown preset <%s>, using default basis: <karl 4>\n", presetstr);
            kmerlen = 4;
            *bf = karsn;
        }
    }
    else {
        basisstr = args_getstr(args, "basis", 0, "karl");
        maxi = sizeof(possible_bases) / sizeof(char*);
        for (i = 0;i < maxi;i++) {
            if (strcmp(possible_bases[i], basisstr) == 0)break;
        }
        if (i >= maxi) {
            if (beginswith("karl", basisstr)) {
                if (basisstr[4] == '*') {
                    kmerlen = atoll(basisstr + 5);
                    *bf = karsn;
                    args_report_warning(NULL, "'-b karl*%d' is deprecated, use '-b karl %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
                else {
                    kmerlen = atoll(basisstr + 4);
                    *bf = karln;
                    args_report_warning(NULL, "'-b karl%d' is deprecated, use '-b s_karl %d' instead\n",(int)kmerlen, (int)kmerlen);
                }
            }
            else if (beginswith("freq", basisstr)) {
                if(basisstr[4] == '*') {
                    kmerlen = atoll(basisstr + 5);
                    *bf = freqs;
                    args_report_warning(NULL, "'-b freq*%d' is deprecated, use '-b freq %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
                else {
                    kmerlen = atoll(basisstr + 4);
                    *bf = freqn;
                    args_report_warning(NULL, "'-b freq%d' is deprecated, use '-b s_freq %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
            }
            else if (beginswith("markz", basisstr)) {
                if (basisstr[5] == '*') {
                    kmerlen = atoll(basisstr + 6);
                    *bf = markz2;
                    args_report_warning(NULL, "'-b markz*%d' is deprecated, use '-b exp %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
                else {
                    kmerlen = atoll(basisstr + 5);
                    *bf = markz1;
                    args_report_warning(NULL, "'-b markz%d' is deprecated, use '-b s_exp %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
            }
            else if (strcmp("TETRA", basisstr) == 0) {
                *bf = TETRA;
                kmerlen = 4;
                args_report_warning(NULL, "'-b TETRA' is deprecated, use '-p TETRA' instead\n");
            }
            else {
                args_report_warning(NULL, "Unknown basis <%s>, defaulting to <karl>\n", basisstr);
                *bf = karsn;
                basisstr = "karl";
            }
        }
        else {
            *bf = bf_list[i];
        }
        if (strcmp(basisstr, "GC") == 0)
            kmerlen = args_getint(args, "basis", 0, 0);
        else
            kmerlen = args_getint(args, "basis", 0, 4);
        if (kmerlen < 0) kmerlen = 0;
        *sig_len = (1LL << (2 * kmerlen));
        args_report_info(NULL, "Calculating length-%d %s signatures\n", (int)kmerlen, basisstr);
    }
    *kmer_len = kmerlen;
    return nocalc;
}
int set_method_function(args_t* args, method_function* mf, double* metharg) {
    char* methstr;
    char* presetstr;
    double tmp;
    presetstr = args_getstr(args, "preset", 0, NULL);
    if (presetstr) {
        if (strcmp(presetstr, "TETRA") == 0) {
            *mf = corr;
            *metharg = 0.0;
        }
        else if (strcmp(presetstr, "PaSiT4") == 0) {
            *mf = SVC;
            *metharg = 0.02;
        }
        else if (strcmp(presetstr, "PaSiT6") == 0) {
            *mf = SVC;
            *metharg = 0.02;
        }
        else {
            args_report_warning(NULL, "Unknown preset <%s>, using default method: <PaSiT 0.02>\n", presetstr);
            *mf = SVC;
            *metharg = 0.02;
        }
    }
    else {
        methstr = args_getstr(args, "method", 0, "SVC");
        if (strcmp(methstr, "AMD") == 0 || strcmp(methstr, "manhattan") == 0) {
            *mf = AMD;
            *metharg = 0.0;
        }
        else if (strcmp(methstr, "ED") == 0 || strcmp(methstr, "euclid") == 0) {
            *mf = ED;
            *metharg = 0.0;
        }
        else if (strcmp(methstr, "SVC") == 0 || strcmp(methstr, "PaSiT") == 0) {
            *mf = SVC;
            *metharg = args_getdouble(args, "method", 0, 0.02);
        }
        else if (strcmp(methstr, "hamming") == 0) {
            *mf = SVC;
            *metharg = 0.0;
        }
        else if (strcmp(methstr, "corr") == 0 || strcmp(methstr, "pearson") == 0) {
            *mf = corr;
            *metharg = 0.0;
        }
        else if (beginswith("SVC", methstr)) {
            *metharg = atof(methstr + 3);
            *mf = SVC;
            args_report_warning(NULL, "'-m SVC<f>' is deprecated, use '-m SVC <f>' instead\n");
        }
        else {
            tmp = 0.0;
            *metharg = args_getdouble(args, "method", 0, tmp);
            if (*metharg != tmp) {
                args_report_warning(NULL, "Unknown method <%s %f>, using default: <PaSiT 0.02>\n", presetstr);
            }
            else {
                args_report_warning(NULL, "Unknown method <%s>, using default: <PaSiT 0.02>\n", presetstr);
            }
            *mf = SVC;
            *metharg = 0.02;
        }
    }
    return 0;
}

double* readsig_fromfile(PF_t* source, size_t* sl) {
    double* data;
    int32_t type;
    int32_t datatype; /* NOTE: datatype is ignored - values are assumed to be stored as double */
    uint64_t n;
    /* it is assumed that the file's endianness matches that of the OS */
    PFread(&type, sizeof(int32_t), 1, source);
    PFread(&datatype, sizeof(int32_t), 1, source);
    PFread(&n, sizeof(uint64_t), 1, source);
    data = malloc(sizeof(double)*n);
    PFread(data, sizeof(double), n, source);
    *sl = (size_t)n;
    return data;
}
signature_t* read_sig_file(char* filename) {
    signature_t* res;
    PF_t* f;
    size_t n;
    res = signature_alloc();
    PFopen(&f, filename, "r");
    if (f) {
        res->fname = os_rmdirname(filename);
        res->data = (void*)readsig_fromfile(f, &n);
        res->len = n;
        PFclose(f);
    }
    return res;
}

signature_t* get_signature_from_fasta(char* filename, basis_function basis, int nlen) {
    PF_t* f;
    nucseq** nsa;
    nucseq tmpseq = EMPTYSEQ;
    nucseq spacerseq = EMPTYSEQ;
    signature_t* sg;
    size_t outcount;
    size_t i;
    PFopen(&f, filename, "r");
    if (!f) {
        args_report_info(NULL, "Failed to read: %s\n", filename);
        return NULL;
    }
    nsa = nucseq_array_from_fasta(f, &outcount, 1, 1);
    nucseq_from_string(&spacerseq, "N");

    sg = signature_alloc();
    sg->fname = os_rmdirname(filename);
    clear_nucseq(&tmpseq);
    for (i = 0;i < outcount;i++) {
        appendtoseq(&tmpseq, &spacerseq);
        appendtoseq(&tmpseq, nsa[i]);
    }
    sg->len = basis(&(tmpseq), nlen, (double**)(&(sg->data)));
    clear_nucseq(&tmpseq);
    clear_nucseq(&spacerseq);
    for (i = 0;i<outcount;i++) {
        clear_nucseq(nsa[i]);
        free(nsa[i]);
    }
    PFclose(f);
    free(nsa);
    return sg;
}
signature_t* get_signatures(char* filename, basis_function basis, int nlen) {
    if (endswith(".sig", filename))
        return read_sig_file(filename);
    else
        return get_signature_from_fasta(filename, basis, nlen);

}

signature_t** GenDisCal_get_signatures(args_t* args, char** sourcefiles, size_t nfiles, size_t* sig_len, int numthreads) {
    basis_function bf;
    size_t i, n;
    size_t* p_n;
    size_t notax;
    size_t* p_notax;
    size_t countsperproc;
    int nocalc;
    size_t kmerlen;
    taxextractor_t* te;
    char* line;
    char* taxtype;
    char* badload;
    char* badtax;
    PF_t* f;
    char** taxstrings;
    size_t ntaxstrings;
    size_t nalloctaxstrings;
    DM64_t* file2tax;
    signature_t** result;
    
    nocalc = set_basis_function(args, &bf, sig_len, &kmerlen);

    result = (signature_t**)calloc(nfiles, sizeof(signature_t*));
    if (!nocalc) {
        /* prepare taxonomy */
        ntaxstrings = 1;
        nalloctaxstrings = 256;
        taxstrings = calloc(nalloctaxstrings, sizeof(char*));
        file2tax = new_DM64(DM_ALGORITHM_BASICSORTEDLIST, 0);
        te = taxextractor_alloc();
        taxtype = NULL;
        if (args_ispresent(args, "taxonomy")) {
            line = args_getstr(args, "taxonomy", 0, "taxonomy.csv");
            taxtype = args_getstr(args, "taxonomy", 1, "PCOF_S");
            if (PFopen(&f, line, "r")) args_report_warning(NULL, "Taxonomy file <%s> could not be loaded\n", line);
            else {
                args_report_progress(NULL, "Loading taxonomy...\n");
                while (line = PFreadline(f)) {
                    if (ntaxstrings >= nalloctaxstrings) {
                        nalloctaxstrings *= 2;
                        taxstrings = (char**) realloc(taxstrings, sizeof(char*)*nalloctaxstrings);
                    }
                    n = text_nextchar(line, 0, ',');
                    taxstrings[ntaxstrings] = (char*)malloc(strlen(line + n + 1) + 1);
                    memcpy(taxstrings[ntaxstrings], line + n + 1, strlen(line + n + 1) + 1);
                    line[n] = 0;
                    taxextractor_add(te, taxtype, line + n + 1);
                    DM64_append(file2tax, line, (int)strlen(line), ntaxstrings);
                    free(line);
                    ntaxstrings++;
                }
                /* te should be sorted before going into the parallel section */
                taxextractor_sort(te); 
                args_report_progress(NULL, "Finished loading taxonomy\n");
            }
        }

        /* read the signatures */
        args_report_progress(NULL, "Reading signatures...\n");
        if (nfiles < (size_t) numthreads) numthreads = (int)nfiles;
        omp_set_num_threads(numthreads);
        n = 0;
        p_n = &n;
        notax = 0;
        p_notax = &notax;
        badload = NULL;
        badtax = NULL;
        countsperproc = (nfiles / numthreads);
        if (countsperproc*numthreads < nfiles)
            countsperproc++;
        args_report_info(NULL, "Each thread handles " _LLD_ " files\n", countsperproc);
        #ifndef NOOPENMP
        #pragma omp parallel
        #endif
        {
            char* ltaxonstring;
            int nullflag;
            size_t ifl_added = 0;
            size_t taxnameid;
            int cthread = omp_get_thread_num();
            size_t firstf = countsperproc*cthread;
            size_t lastf = countsperproc*(cthread + 1);
            if (firstf > nfiles) firstf = nfiles;
            if (lastf > nfiles) lastf = nfiles;
            #pragma omp critical
            args_report_info(NULL, "Thread %3d handles [%zd : %zd]\n", cthread, firstf, lastf);
            for (ifl_added = firstf;ifl_added < lastf;ifl_added++) {
                result[ifl_added] = get_signatures(sourcefiles[ifl_added], bf, (int)kmerlen);
                if (taxtype) {
                    nullflag = 0;
                    taxnameid = DM64_get(file2tax, result[ifl_added]->fname, (int)strlen(result[ifl_added]->fname), &nullflag);
                    ltaxonstring = taxstrings[taxnameid];
                    if (!nullflag)
                        taxextractor_translate(te, taxtype, ltaxonstring, result[ifl_added]->lineage, 9);
                    else {
                        #ifndef NOOPENMP
                        #pragma omp atomic
                        #endif
                        (*p_notax)++;
                        badtax = sourcefiles[ifl_added];
                    }
                }
                if (result[ifl_added]) {
                    #ifndef NOOPENMP
                    #pragma omp atomic
                    #endif
                    (*p_n)++;
                }
                else {
                    badload = sourcefiles[ifl_added];
                }
                if (nfiles < 100 || (*p_n) % (nfiles / 100) == 0) {
                    #pragma omp critical
                    args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files loaded\r", n, nfiles, (int)((n * 100) / nfiles));
                }
            }
        }
        #ifndef NOOPENMP
        omp_set_num_threads(1);
        #endif
        args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files loaded \n", n, nfiles, (int)((n * 100) / nfiles));
        if (n != nfiles) {
            if (nfiles - n != 1)
                args_report_warning(NULL, "Some " _LLD_ " files could not be properly loaded\n", nfiles - n);
            else if(badload)
                args_report_warning(NULL, "<%> could not be properly loaded\n",badload);
        }
        if (notax > 0) {
            if (notax != 1)
                args_report_warning(NULL, "Taxonomy could not be assigned for " _LLD_ " files\n", nfiles - n);
            else if(badtax)
                args_report_warning(NULL, "Taxonomy could not be assigned for <%>\n",badtax);
        }
        args_report_progress(NULL, "Finished reading signatures\n");
        /* cleanup */
        free_DM64(file2tax);
        for (i = 0;i < ntaxstrings;i++) {
            if (taxstrings[i])free(taxstrings[i]);
        }
        free(taxstrings);
        taxextractor_free(te);
        args_report_info(NULL, "Signature cleanup complete.\n");
    }
    
    return result;
}

#define OUTPUTMODE_NORMAL   0
#define OUTPUTMODE_MATRIX   1
#define OUTPUTMODE_HIST     2
int strtotaxsim(const char* str) {
    if (strcmp(str, "?") == 0)return 0;
    if (strcmp(str, "Different_Phyla") == 0)return 0;
    if (strcmp(str, "Same_Phylum") == 0)return 1;
    if (strcmp(str, "Same_Class") == 0)return 2;
    if (strcmp(str, "Same_Order") == 0)return 3;
    if (strcmp(str, "Same_Family") == 0)return 4;
    if (strcmp(str, "Same_Genus") == 0)return 5;
    if (strcmp(str, "Same_Species") == 0)return 6;
    if (strcmp(str, "Same_Subspecies") == 0)return 7;
    if (strcmp(str, "Replicate") == 0)return 8;
    return -1;
}
char* taxsimtostr(int val) {
    switch (val) {
    case 0:
        return "Different_Phyla";
    case 1:
        return "Same_Phylum";
    case 2:
        return "Same_Class";
    case 3:
        return "Same_Order";
    case 4:
        return "Same_Family";
    case 5:
        return "Same_Genus";
    case 6:
        return "Same_Species";
    case 7:
        return "Same_Subspecies";
    case 8:
        return "Replicate";
    default:
        return "?";
    }
}


int GenDisCal_perform_comparisons(args_t* args, signature_t** signatures, size_t nfiles, size_t sig_len, int numthreads) {
    method_function mf;
    double metharg;
    int nocalc;
    int outputmode;
    double binwidth;
    double minval;
    double maxval;
    int onlylevel;
    int issearch;
    datatable_t* otable;
    PF_t* outputf;
    size_t i,n;
    size_t* p_n;
    size_t nfsquare;
    size_t countsperproc;
    
    /* parse relevant arguments */
    nocalc = set_method_function(args, &mf, &metharg);
    args_report_info(NULL, "Method function selected successfully\n", numthreads);
    
    outputmode = OUTPUTMODE_NORMAL;
    if (args_ispresent(args, "distancematrix"))
        outputmode = OUTPUTMODE_MATRIX;
    if (args_ispresent(args, "histogram")) {
        outputmode = OUTPUTMODE_HIST;
        binwidth = args_getdouble(args, "histogram", 0, 0.001);
    }

    if (outputmode == OUTPUTMODE_NORMAL) {
        args_report_info(NULL, "Output mode: NORMAL\n");
    }
    if (outputmode == OUTPUTMODE_HIST) {
        args_report_info(NULL, "Output mode: HISTOGRAM\n");
    }
    if (outputmode == OUTPUTMODE_MATRIX) {
        args_report_info(NULL, "Output mode: DISTANCE MATRIX\n");
    }

    minval = args_getdouble(args, "above", 0, -INFINITY);
    maxval = args_getdouble(args, "below", 0, INFINITY);
    if (args_ispresent(args, "only"))
        onlylevel = strtotaxsim(args_getstr(args, "only", 0, "Same_Species"));
    else
        onlylevel = -2;
    args_report_info(NULL, "Limiting results to: [%f:%f] and only %d\n", minval, maxval, onlylevel);
    if (args_ispresent(args, "search")) {
        issearch = 1;
        if (outputmode == OUTPUTMODE_MATRIX)
            outputmode = OUTPUTMODE_NORMAL;
        nfiles--;
        args_report_info(NULL, "Search mode is enabled\n");
    }
    else issearch = 0;

    PFopen(&outputf, args_getstr(args, "output", 0, "stdout"),"w");
    if (!outputf) {
        args_report_error(NULL, "The output file <%s> could not be opened for writing! Aborting.\n", args_getstr(args, "output", 0, "stdout"));
        return 1;
    }

    /* create the output table if needed */

    otable = NULL;
    if (outputmode == OUTPUTMODE_MATRIX) {
        otable = datatable_alloc(nfiles, nfiles, NAN);
        for (i = 0;i < nfiles;i++) {
            if (signatures[i] && signatures[i]->fname) {
                datatable_setcolname(otable, i, signatures[i]->fname);
                datatable_setrowname(otable, i, signatures[i]->fname);
            }
            else {
                datatable_setcolname(otable, i, "NA");
                datatable_setrowname(otable, i, "NA");
            }
        }
    }
    if (outputmode == OUTPUTMODE_HIST) {
        /*otable = datatable_alloc(1+10, 1+(size_t)(1.0/binwidth), 0.0);*/
        otable = datatable_alloc(1, 1, 0.0);
        datatable_setcolname(otable, 0, "bin");
        datatable_set(otable, 0, 0, binwidth / 2.0);
    }
    args_report_info(NULL, "Argument parsing is done. Beginning comparisons.\n");

    /* do the actual calculations */
    if (!nocalc) {
        if (nfiles < (size_t) numthreads) numthreads = (int) nfiles;
        omp_set_num_threads(numthreads);
        n = 0;
        p_n = &n;
        nfsquare = nfiles*nfiles/2;
        countsperproc = (nfiles / numthreads);
        if (countsperproc*numthreads < nfiles)
            countsperproc++;
        args_report_info(NULL, "Each thread handles " _LLD_ " files\n", countsperproc);
        #ifndef NOOPENMP
        #pragma omp parallel
        #endif
        {
            char* fn1;
            char* fn2;
            double dist = 0.0;
            double oldval;
            int simlevel;
            int nullflag;
            size_t ifl_added = 0;
            size_t ifl_c;
            int cthread = omp_get_thread_num();
            size_t firstf = cthread;
            size_t lastf = nfiles;
            size_t tmp1;
            if (firstf > nfiles) firstf = nfiles;
            for (ifl_added = firstf;ifl_added < lastf;ifl_added += numthreads) {
                if (outputmode == OUTPUTMODE_MATRIX) {
                    datatable_set(otable, ifl_added, ifl_added, 0.0);
                }
                if (issearch) {
                    if (signatures[ifl_added] && signatures[nfiles])
                        dist = mf((double*)(signatures[ifl_added]->data), (double*)(signatures[nfiles]->data), sig_len, metharg);
                    else
                        dist = 1.0;
                    simlevel = signature_taxsimlevel(signatures[ifl_added], signatures[nfiles]);
                    if ((onlylevel < -1 || simlevel == onlylevel) && dist >= minval && dist <= maxval) {
                        if (outputmode == OUTPUTMODE_NORMAL) {
                            if (signatures[ifl_added] && signatures[ifl_added]->fname) fn1 = signatures[ifl_added]->fname;
                            else fn1 = "<File could not opened>";
                            if (signatures[nfiles] && signatures[nfiles]->fname) fn2 = signatures[nfiles]->fname;
                            else fn2 = "<File could not opened>";
                            #pragma omp critical
                            PFprintf(outputf, "%s,%s,%s,%f\n", fn1, fn2, taxsimtostr(simlevel), dist);
                        }
                        else if (outputmode == OUTPUTMODE_HIST) {
                            if (dist >= 0.0) {
                                tmp1 = (size_t)round(dist / binwidth);
                                #pragma omp critical
                                {
                                    nullflag = 0;
                                    oldval = datatable_get(otable, simlevel + 2, tmp1, &nullflag);
                                    if (nullflag)
                                        oldval = 0.0;
                                    datatable_set(otable, simlevel + 2, tmp1, oldval + 1.0);
                                }
                            }
                        }
                        else if (outputmode == OUTPUTMODE_MATRIX) {
                            datatable_set(otable, ifl_added, nfiles, dist);
                            datatable_set(otable, nfiles, ifl_added, dist);
                        }
                        #pragma omp atomic
                        (*p_n)++;
                        if (nfsquare < 100 || (*p_n) % (nfiles / 100) == 0) {
                            #pragma omp critical
                            args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files compared\r", (*p_n), nfiles, (int)(((double)n * 100) / (double)nfiles));
                        }
                    }
                }
                else {
                    for (ifl_c = ifl_added + 1; ifl_c < nfiles;ifl_c++) {
                        if (signatures[ifl_added] && signatures[ifl_c])
                            dist = mf((double*)(signatures[ifl_added]->data), (double*)(signatures[ifl_c]->data), sig_len, metharg);
                        else
                            dist = 1.0;
                        simlevel = signature_taxsimlevel(signatures[ifl_added], signatures[ifl_c]);
                        if ((onlylevel < -1 || simlevel == onlylevel) && dist >= minval && dist <= maxval) {
                            if (outputmode == OUTPUTMODE_NORMAL) {
                                if (signatures[ifl_added] && signatures[ifl_added]->fname) fn1 = signatures[ifl_added]->fname;
                                else fn1 = "<File could not opened>";
                                if (signatures[ifl_c] && signatures[ifl_c]->fname) fn2 = signatures[ifl_c]->fname;
                                else fn2 = "<File could not opened>";
                                #pragma omp critical
                                PFprintf(outputf, "%s,%s,%s,%f\n", fn1, fn2, taxsimtostr(simlevel), dist);
                            }
                            else if (outputmode == OUTPUTMODE_HIST) {
                                if (dist >= 0.0) {
                                    tmp1 = (size_t)floor(dist / binwidth);
                                    #pragma omp critical
                                    {
                                        nullflag = 0;
                                        oldval = datatable_get(otable, simlevel + 2, tmp1, &nullflag);
                                        if (nullflag)
                                            oldval = 0.0;
                                        datatable_set(otable, simlevel + 2, tmp1, oldval + 1.0);
                                    }
                                }
                            }
                            else if (outputmode == OUTPUTMODE_MATRIX) {
                                #pragma omp critical
                                {
                                    datatable_set(otable, ifl_added, ifl_c, dist);
                                    datatable_set(otable, ifl_c, ifl_added, dist);
                                }
                            }
                            #pragma omp atomic
                            (*p_n)++;
                            if (nfsquare < 100 || (*p_n) % (nfsquare / 100) == 0) {
                                #pragma omp critical
                                args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) comparisons performed\r", (*p_n), nfsquare, (int)(((double)n * 100) / (double)nfsquare));
                            }
                        }
                    }
                }
            }
        }
        args_report_progress(NULL, "A Total of " _LLD_ " comparisons were performed                        \n", n);
    }
    if (otable) {
        if (outputmode == OUTPUTMODE_HIST) {
            n = datatable_ncols(otable) - 1;
            for (i = 0;i < n;i++) {
                datatable_setcolname(otable, i + 1, taxsimtostr((int)i - 1));
            }
            n = datatable_nrows(otable);
            for (i = 0;i < n;i++) {
                datatable_set(otable, 0, i, binwidth*((double)i+0.5));
            }
            datatable_write(otable, outputf, ",", DATATABLE_WRITECOLNAMES);
        }
        if (outputmode == OUTPUTMODE_MATRIX) {
            datatable_write(otable, outputf, ",", DATATABLE_WRITECOLNAMES | DATATABLE_WRITEROWNAMES);
        }
    }
    /* cleanup */
    PFclose(outputf);
    if (otable)datatable_free(otable);

    return 0;
}

int GenDisCal(args_t* args) {
    char** filelist;
    signature_t** siglist;
    size_t i;
    size_t numfiles;
    size_t siglen;
    int numthreads;
    filelist = GenDisCal_get_file_list(args, &numfiles);
    if (numfiles < 2) {
        args_report_error(NULL, "Number of files (" _LLD_ ") to be loaded is below 2! Aborting.\n", numfiles);
    }
    else {
        numthreads = args_getint(args, "nprocs", 0, autothreads());
        args_report_info(NULL, "Using up to %d threads\n", numthreads);
        siglist = GenDisCal_get_signatures(args, filelist, numfiles, &siglen, numthreads);
        GenDisCal_perform_comparisons(args, siglist, numfiles, siglen, numthreads);

        for (i = 0;i < numfiles;i++) {
            if (siglist[i])signature_free(siglist[i]);
            siglist[i] = NULL;
        }
        free(siglist);
    }
    for (i = 0;i < numfiles;i++) {
        free(filelist[i]);
    }
    free(filelist);
    return 0;
}

#ifdef TEST
#else
int main(int argc, char** argv) {
    args_t* args;
    int result;
    args = GenDisCal_init_args(argc, argv);
    if(!args_is_helpmode(args))
        result = GenDisCal(args);
    args_free(args);
    return result;
}
#endif
