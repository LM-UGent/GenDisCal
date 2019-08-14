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
    args_add(result, "ordered", 'd', "");
    args_add(result, "sort_by_distance", 'r', "");
    args_add(result, "below", 'w', "float");
    args_add(result, "above", 'v', "float");
    args_add(result, "only", 'y', "str");
    args_add(result, "filtershort", 'l', "int");
    args_add(result, "keepsigs", 'k',"str");
#ifndef NOOPENMP
    args_add(result, "nprocs", 'n', "int");
#endif

    args_add_help(result, "basis", "SIGNATURE BASIS",
        "Signature bases are used to transform a genome into a vector of numbers. "
        "Available signatures are as follows:\n"
        "   len        -- use genome length as a signature\n"
        "   GC         -- use GC content as a signature\n"
        "   freq <k>   -- use frequencies of <k>-mers as a signature\n"
        "   karl <k>   -- use karlin signatures of <k>-mers\n"
        "   exp <k>    -- use <k>-mer bias with regards to monomer composition\n"
        /*"   bias <k>   -- use <k>-mer bias with regards to <k-1>-mer composition\n"*/
        "For single-stranded DNA and RNA, alternative versions of these functions are provided, "
        "beginning with 's_', as follows:\n"
        "   s_freq, s_karl, s_exp\n",
        "signature basis to be used");
    args_add_help(result, "method", "DISTANCE COMPUTATION METHOD",
        "Distance computation methods are used to compute the dissimilarity between genome signatures."
        "Available methods are as follows:\n"
        "   manhattan/AMD  -- use the manhattan distance\n"
        "   euclid/ED      -- use the euclidian distance\n"
        "   hamming        -- use the hamming distance\n"
        "   PaSiT/SVC <f>  -- use the PaSiT distance with a non-standard threshold\n"
        "   pearson/corr   -- use Pearson correlation as dissimilarity measure\n"
        "   reldist        -- use the mean of pairwise relative distances between signature values\n",
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
        " > approxANI (-b bias 4 -m corr)\n"
        "    approxANI is a minhash-based approximation of Average Nucleotide\n"
        "    Identity (ANI). This method is similar to the one used by the \n"
        "    Mash software (Ondov et al., Gen. biol., 2016) "
        "    The ambiguity region for species is [0.04-0.06]\n"
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
        /*"  semicol -- filename;Phylum;Class;Order;Family;Genus;species;etc\n"*/
        /*"  pipe    -- filename | Phylum | Class | Order | Family | Genus | Species | etc\n"*/,
        "taxonomy file");
    args_add_help(result, "ordered", "FORCE ORDERED OUTPUT",
        "If this option is present, the output will be forcefully ordered. This has no effect "
        "on histogram and distance matrix output, but will ensure that the order in normal "
        "output mode stays consistent across multiple runs when using more than one thread.\n"
        "This comes at the cost of higher RAM usage.\n",
        "force ordered output");
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
    args_add_help(result, "sort_by_distance", "SORT BY DISTANCE",
        "Sort output table according to distance. Automatically sets the --ordered flag when used. "
        "Ignored when --distancematrix or --histogram is specified (default: no)\n",
        "Sort output table according to distance");
    args_add_help(result, "below", "VALUES BELOW A GIVEN VALUE",
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
    args_add_help(result, "filtershort", "FILTER SHORT SEQUENCES",
        "Indicates the minimum contig length to keep. Experience shows that contigs which "
        "are shorther than ~2000 bp vary depending on the assembly method. (default: 2000)",
        "Indicates the minimum contig length to keep.");
    args_add_help(result, "keepsigs", "KEEP SIGNATURES",
        "If this option is present, a signature file will be generated for re-use, named "
        "<filename>.sig for individual signature files, and <filename>.sdb for multiple "
        "signature files. The argument specifies whether or not the comparisons should "
        "still be performed after this step:\n"
        "  continue : perform comparisons after creating signature files (default)\n"
        "  stop     : abort the program creating signature files\n",
        "generate signature files");
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
    size_t multicount;
    int datatype;
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
    size_t i;
    signature_t** demux;
    if (target) {
        if (target->multicount > 0) {
            demux = (signature_t**)(target->data);
            for (i = 0;i < target->multicount;i++) {
                signature_free(demux[i]);
            }
            target->multicount = 0;
        }
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
signature_t* signature_multiplex(signature_t** signatures, size_t count) {
    signature_t* res;
    res = signature_alloc();
    res->data = signatures;
    res->multicount = count;
    return res;
}
signature_t** signature_demultiplex(signature_t* source) {
    signature_t** res;
    res = (signature_t**)(source->data);
    source->data = NULL;
    source->multicount = 0;
    source->len = 0;
    return res;
}
int signature_is_multiplexed(signature_t* target) {
    return target->multicount != 0;
}
size_t signature_get_multiplex_amount(signature_t* target) {
    return target->multicount;
}
int set_basis_function(args_t* args, basis_function* bf, size_t* sig_len, size_t* kmer_len) {
    size_t maxi, i;
    size_t kmerlen;
    int nocalc;
    char* basisstr;
    char* presetstr;
    char* possible_bases[] = { "GC","freq","karl","exp","s_freq","s_karl","s_exp","len","SSU" };
    basis_function bf_list[] = { freqs,freqs,karsn,markz2,freqn,karln,markz1,gensz,SSUseq };
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
        else if (strcmp(presetstr, "approxANI") == 0) {
            *bf = minhashsig;
            kmerlen = 21;
            *sig_len = SIGLEN_MINHASH;
        }
        else if (strcmp(presetstr, "combinedSpecies") == 0) {
            *bf = combinedsig;
            kmerlen = 21;
            *sig_len = SIGLEN_MINHASH;
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
            kmerlen = 0;
        else if (strcmp(basisstr, "len") == 0)
            kmerlen = 0;
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
        else if (strcmp(presetstr, "approxANI") == 0) {
            *mf = approxANI;
            *metharg = 21.0;
        }
        else if (strcmp(presetstr, "combinedSpecies") == 0) {
            *mf = combinedSpecies;
            *metharg = 21.0;
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
        else if (strcmp(methstr, "reldist") == 0) {
            *mf = reldist;
            *metharg = 0.0;
        }
        else if (beginswith("SVC", methstr)) {
            *metharg = atof(methstr + 3);
            *mf = SVC;
            args_report_warning(NULL, "'-m SVC<f>' is deprecated, use '-m SVC <f>' instead\n");
        }
        else if (beginswith("localalign", methstr)) {
            *metharg = 1;
            *mf = localalign;
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

void read_sigmeta_from_file(PF_t* source, signature_t* sig) {
    int64_t fnamelen;
    /* it is assumed that the file's endianness matches that of the OS */
    fnamelen = PFgetint64(source);
    if (fnamelen > 0) {
        sig->fname = malloc(sizeof(char)*(fnamelen + 1));
        sig->fname[fnamelen] = 0;
        sig->freefname = 1;
        PFread(sig->fname, sizeof(char), fnamelen, source);
    }
}
signature_t* readsig_fromfile(PF_t* source) {
    signature_t* res;
    double* data;
    int32_t type;
    int32_t datatype;
    uint64_t n;
    /* it is assumed that the file's endianness matches that of the OS */
    PFread(&type, sizeof(int32_t), 1, source);
    PFread(&datatype, sizeof(int32_t), 1, source);
    PFread(&n, sizeof(uint64_t), 1, source);
    if (datatype == 2) {
        /* rem: datatype == 2 is bytes */
        data = malloc(n);
        PFread(data, 1, n, source);
    }
    else {
        data = malloc(sizeof(double)*n);
        PFread(data, sizeof(double), n, source);
    }
    res = signature_alloc();
    res->data = data;
    res->len = n;
    res->datatype = datatype;
    read_sigmeta_from_file(source, res);
    return res;
}
signature_t* read_sig_file(char* filename) {
    signature_t* res;
    PF_t* f;
    PFopen(&f, filename, "rb");
    if (f) {
        res = readsig_fromfile(f);
        if (!res->fname)
            res->fname = os_rmdirname(filename);
        PFclose(f);
    }
    else
        res = signature_alloc();
    return res;
}
signature_t* read_multisig_file_as_multiplexed(char* filename, int64_t maxcount) {
    signature_t* res;
    signature_t** resarray;
    PF_t* f;
    int64_t cursig;
    int64_t fcount;
    res = NULL;
    PFopen(&f, filename, "rb");
    if (f) {
        fcount = (size_t)PFgetint64(f);
        resarray = calloc((size_t)fcount, sizeof(signature_t*));
        if (fcount > maxcount)fcount = maxcount;
        for (cursig = 0; cursig < fcount;cursig++) {
            res = readsig_fromfile(f);
            resarray[cursig] = res;
        }
        res = signature_multiplex(resarray, fcount);
        PFclose(f);
    }
    return res;
}
signature_t* get_signature_from_fasta(char* filename, basis_function basis, int nlen, size_t minseqlen) {
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
    nsa = nucseq_array_from_fasta(f, &outcount, 1, minseqlen);
    nucseq_from_string(&spacerseq, "N");

    sg = signature_alloc();
    sg->fname = os_rmdirname(filename);
    clear_nucseq(&tmpseq);
    for (i = 0;i < outcount;i++) {
        appendtoseq(&tmpseq, &spacerseq);
        appendtoseq(&tmpseq, nsa[i]);
    }
    sg->len = basis(&(tmpseq), nlen, (double**)(&(sg->data)));
    sg->datatype = 2; /* starting from version 0.2, all computed signatures are to be considered as binary data */
    clear_nucseq(&tmpseq);
    clear_nucseq(&spacerseq);
    for (i = 0;i<outcount;i++) {
        clear_nucseq(nsa[i]);
        free(nsa[i]);
    }
    PFclose(f);
    free(nsa);
    if (!sg->data) {
        signature_free(sg);
        sg = NULL;
    }
    return sg;
}
signature_t* get_signatures(char* filename, basis_function basis, int nlen, size_t minseqlen) {
    if (basis == SSUseq) {
        if (endswith(".sig", filename) || endswith(".sdb", filename)) {
            args_report_warning(NULL, "SSUs cannot be stored as signatures. Skipping file <%s>\n", filename);
            return NULL;
        }
        else {
            return get_signature_from_fasta(filename, basis, nlen, minseqlen);
        }
    }
    else {
        if (endswith(".sig", filename))
            return read_sig_file(filename);
        else if (endswith(".sdb", filename)) {
            return read_multisig_file_as_multiplexed(filename, 100000000);
        }
        else
            return get_signature_from_fasta(filename, basis, nlen, minseqlen);

    }
}

void write_signature(PF_t* target, signature_t* source) {
    int32_t type;
    int32_t datatype; /* DATATYPE: 1 - double, 2 - byte */
    uint64_t fnamelen;
    datatype = source->datatype;
    type = 1;
    /* it is assumed that the file's endianness matches that of the OS */
    PFwrite(&type, sizeof(int32_t), 1, target);
    PFwrite(&datatype, sizeof(int32_t), 1, target);
    PFwrite(&(source->len), sizeof(uint64_t), 1, target);
    if (datatype <= 1)
        PFwrite(source->data, sizeof(double), source->len, target);
    else if (datatype == 2)
        PFwrite(source->data, 1, source->len, target);
    if (source->fname) {
        /* inserted a the end to maintain portability */
        fnamelen = strlen(source->fname);
        PFputint64(target,fnamelen);
        PFwrite(source->fname, 1, strlen(source->fname), target);
    }
}
void write_sigfile(const char* basename, signature_t* source) {
    PF_t* f;
    char* newname;
    size_t baselen;
    baselen = strlen(basename);
    newname = (char*) malloc(baselen + 5);
    memcpy(newname, basename, baselen);
    newname[baselen] = '.';
    newname[baselen + 1] = 's';
    newname[baselen + 2] = 'i';
    newname[baselen + 3] = 'g';
    newname[baselen + 4] = '\0';
    PFopen(&f, newname, "wb");
    write_signature(f, source);
    PFclose(f);
}
void write_multisig_file(const char* basename, signature_t** signatures, size_t nsigs) {
    PF_t* f;
    char* newname;
    size_t baselen;
    size_t n;
    baselen = strlen(basename);
    newname = (char*)malloc(baselen + 5);
    memcpy(newname, basename, baselen);
    newname[baselen] = '.';
    newname[baselen + 1] = 's';
    newname[baselen + 2] = 'd';
    newname[baselen + 3] = 'b';
    newname[baselen + 4] = '\0';
    PFopen(&f, newname, "wb");
    PFputint64(f, (int64_t)nsigs);
    for (n = 0;n < nsigs;n++) {
        write_signature(f, signatures[n]);
    }
    PFclose(f);
}
signature_t** GenDisCal_get_signatures(args_t* args, char** sourcefiles, size_t* numfiles, size_t* sig_len, int numthreads) {
    basis_function bf;
    size_t i, j, k, n;
    size_t* p_n;
    size_t* i_p;
    size_t notax;
    size_t* p_notax;
    size_t minseqlen;
    size_t countsperproc;
    size_t extrafiles;
    int nocalc;
    int nflag;
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
    signature_t** tmpresult;
    signature_t** subres;
    int64_t taxnameid_demux;
    int gensigs;
    size_t nfiles;
    nfiles = *numfiles;
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
        gensigs = (args_ispresent(args, "keepsigs")?1:0);
        if (gensigs && args_ispresent(args, "filelist"))gensigs = 2;
        args_report_progress(NULL, "Reading signatures...\n");
        if (bf != SSUseq) minseqlen = 2000;
        else minseqlen = 1200;
        minseqlen = (size_t)args_getint(args, "filtershort", 0, (int)minseqlen);
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
        args_report_info(NULL, "Each thread handles on average " _LLD_ " files\n", countsperproc);
        i = 0;
        i_p = &i;
        #pragma omp parallel
        {
            char* ltaxonstring;
            int nullflag;
            int muxed;
            size_t ifl_added = 0;
            size_t taxnameid;
            int cthread = omp_get_thread_num();
            #pragma omp critical
            ifl_added = ((*i_p)++);
            while(ifl_added<nfiles) {
                result[ifl_added] = get_signatures(sourcefiles[ifl_added], bf, (int)kmerlen, minseqlen);
                muxed = signature_is_multiplexed(result[ifl_added]);
                if (!muxed && gensigs==1 && !endswith(".sig", sourcefiles[ifl_added]))
                    write_sigfile(sourcefiles[ifl_added], result[ifl_added]);
                if (!muxed && taxtype) {
                    nullflag = 0;
                    taxnameid = DM64_get(file2tax, result[ifl_added]->fname, (int)strlen(result[ifl_added]->fname), &nullflag);
                    ltaxonstring = taxstrings[taxnameid];
                    if (!nullflag)
                        taxextractor_translate(te, taxtype, ltaxonstring, result[ifl_added]->lineage, 9);
                    else {
                        #pragma omp atomic
                        (*p_notax)++;
                        badtax = sourcefiles[ifl_added];
                    }
                }
                if (result[ifl_added]) {
                    #pragma omp atomic
                    (*p_n)++;
                }
                else {
                    badload = sourcefiles[ifl_added];
                }
                if (nfiles < 100 || (*p_n) % (nfiles / 100) == 0) {
                    #pragma omp critical
                    args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files loaded\r", n, nfiles, (int)((n * 100) / nfiles));
                }
                #pragma omp critical
                ifl_added = ((*i_p)++);
            }
        }
        #ifndef NOOPENMP
        omp_set_num_threads(1);
        #endif
        /* report loading result */
        args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files loaded \n", n, nfiles, (int)((n * 100) / nfiles));
        if (n != nfiles) {
            if (nfiles - n != 1)
                args_report_warning(NULL, "Some " _LLD_ " files could not be properly loaded\n", nfiles - n);
            else if (badload)
                args_report_warning(NULL, "<%> could not be properly loaded\n", badload);
        }
        /* demultiplex multiplexed signatures */
        nfiles = n;
        j = 0;
        for (i = 0;i<n;i++) {
            if (signature_is_multiplexed(result[i])) {
                nfiles += signature_get_multiplex_amount(result[i]) - 1;
                j++;
            }
        }
        if (n != nfiles) {
            args_report_progress(NULL, "Additionally, " _LLD_ " signature(s) are contained inside of " _LLD_ " SDB file(s)\n", nfiles - n + j, j);
            tmpresult = (signature_t**)malloc(sizeof(signature_t*)*nfiles);
            extrafiles = 0;
            for (i = 0;i < n;i++) {
                if (signature_is_multiplexed(result[i])) {
                    j = signature_get_multiplex_amount(result[i]);
                    subres = signature_demultiplex(result[i]);
                    signature_free(result[i]);
                    memcpy(tmpresult + extrafiles, subres, sizeof(signature_t*)*(j));
                    free(subres);
                    if (taxtype) {
                        /* assign taxonomy */
                        for (k = 0;k < j;k++) {
                            nflag = 0;
                            taxnameid_demux = DM64_get(file2tax, tmpresult[extrafiles + k]->fname, (int)strlen(tmpresult[extrafiles + k]->fname), &nflag);
                            line = taxstrings[taxnameid_demux];
                            if (!nflag)
                                taxextractor_translate(te, taxtype, line, tmpresult[extrafiles + k]->lineage, 9);
                            else {
                                notax++;
                                badtax = sourcefiles[i];
                            }
                        }
                    }
                    extrafiles += j;
                }
                else {
                    tmpresult[extrafiles] = result[i];
                    extrafiles++;
                }
            }
            free(result);
            result = tmpresult;
            nfiles = extrafiles;
        }
        /* generate signature file */
        if (gensigs) {
            if (bf == SSUseq) {
                args_report_warning(NULL, "SSU's cannot be exported as signatures - no signature database file generated\n");
            }
            else {
                line = args_getstr(args, "filelist", 0, "all_sigs.list");
                if (strcmp(args_getstr(args, "keepsigs", 0, "continue"), "stop") == 0)
                    line = args_getstr(args, "output", 0, "all_sigs.list");
                write_multisig_file(line, result, nfiles);
            }
        }
        /* report on taxonomy assignment */
        if (notax > 0) {
            if (notax != 1)
                args_report_warning(NULL, "Taxonomy could not be assigned for " _LLD_ " signatures\n", nfiles - n);
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
    *numfiles = nfiles;
    return result;
}

#define OUTPUTMODE_NORMAL   0
#define OUTPUTMODE_MATRIX   1
#define OUTPUTMODE_HIST     2
#define OUTPUTMODE_ORDERED  3
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


char* GenDisCal_ordered_output_gen_line(signature_t** signatures,size_t fileA, size_t fileB, int simlevel) {
    char* filename1;
    char* filename2;
    char* output;
    char* simlevelstr;
    size_t k;
    simlevelstr = taxsimtostr(simlevel);
    if (signatures[fileA] && signatures[fileA]->fname) filename1 = signatures[fileA]->fname;
    else filename1 = "<File could not opened>";
    if (signatures[fileB] && signatures[fileB]->fname) filename2 = signatures[fileB]->fname;
    else filename2 = "<File could not opened>";
    k = strlen(filename1) + strlen(filename2) + strlen(taxsimtostr(simlevel)) + 3;
    output = malloc(k);
#ifdef _WIN32
    sprintf_s(output, k, "%s,%s,%s", filename1, filename2, simlevelstr);
#else
    sprintf(output, "%s,%s,%s", filename1, filename2, simlevelstr);
#endif
    output[k - 1] = 0;
    return output;
}

/* the amount of arguments for this function is a bit excessive*/
void GenDisCal_print_ordered_output(
        size_t orderedamount, int* ocmp_types, double* ocmp_results, size_t nfiles,
        int issearch, int sort_by_distance, PF_t* outputf, signature_t** signatures
    ) {
    size_t i, j, k, n;
    char** outputlines;
    size_t offset;
    PFprintf(outputf, "File1,File2,Expected_Relation,Distance\n");
    n = 0;
    for (j = 0;j < orderedamount;j++) {
        if (ocmp_types[j] != -3)n++;
    }
    outputlines = (char**)calloc(n, sizeof(char*));
    k = 0;
    if (issearch) {
        for (j = 0;j < orderedamount;j++) {
            if (ocmp_types[j] != -3) {
                outputlines[k] = GenDisCal_ordered_output_gen_line(signatures, nfiles, j, ocmp_types[j]);
                ocmp_types[k] = ocmp_types[j];
                ocmp_results[k] = ocmp_results[j];
                k++;
            }
        }
    }
    else {
        for (i = 0;i < nfiles;i++) {
            offset = i*nfiles + 1 - (i + 1)*(i + 2) / 2;
            for (j = i + 1; j < nfiles;j++) {
                if (ocmp_types[offset + j - 1] != -3) {
                    outputlines[k] = GenDisCal_ordered_output_gen_line(signatures, i, j, ocmp_types[offset + j - 1]);
                    ocmp_types[k] = ocmp_types[offset + j - 1];
                    ocmp_results[k] = ocmp_results[offset + j - 1];
                    k++;
                }
            }
        }
    }
    if (sort_by_distance)
        vec_sort_by(outputlines, sizeof(char*), ocmp_results, n);
    for (j = 0;j < n;j++) {
        PFprintf(outputf, "%s,%f\n", outputlines[j], ocmp_results[j]);
        free(outputlines[j]);
    }
    free(outputlines);
    free(ocmp_results);
    free(ocmp_types);
}

datatable_t* GenDisCal_perform_comparisons_hist_thread(
        args_t* args, signature_t** signatures, size_t nfiles, size_t sig_len,
        int threadid, int numthreads, size_t* p_n, method_function mf, double metharg) {
    double dist = 0.0;
    double oldval;
    int simlevel;
    int nullflag;
    size_t ifl_added = 0;
    size_t ifl_c;
    int cthread = threadid;
    size_t firstf = cthread;
    size_t lastf = nfiles;
    size_t tmp1;
    size_t ordered_offset = 0;
    datatable_t* otable;
    size_t nextprint;
    size_t n;
    size_t nfsquare;

    double minval;
    double maxval;
    double binwidth;
    int onlylevel;
    int issearch;

    binwidth = args_getdouble(args, "histogram", 0, 0.001);
    minval = args_getdouble(args, "above", 0, -INFINITY);
    maxval = args_getdouble(args, "below", 0, INFINITY);
    if (args_ispresent(args, "only"))
        onlylevel = strtotaxsim(args_getstr(args, "only", 0, "Same_Species"));
    else
        onlylevel = -2;
    if (args_ispresent(args, "search")) {
        nfiles--;
        issearch = 1;
    }
    else {
        issearch = 0;
    }
    otable = datatable_alloc(1, 1, 0.0);

    if (firstf > nfiles) firstf = nfiles;
    nfsquare = nfiles*(nfiles - 1) / 2;
    nextprint = 0;
    for (ifl_added = firstf;ifl_added < lastf;ifl_added += numthreads) {
        if (issearch) {
            simlevel = signature_taxsimlevel(signatures[ifl_added], signatures[nfiles]);
            if (onlylevel < -1 || simlevel == onlylevel) {
                if (signatures[ifl_added] && signatures[nfiles])
                    dist = mf((double*)(signatures[ifl_added]->data), (double*)(signatures[nfiles]->data), sig_len, metharg);
                else
                    dist = 1.0;
                if (dist >= minval && dist <= maxval) {
                    if (dist >= 0.0) {
                        tmp1 = (size_t)round(dist / binwidth);
                        nullflag = 0;
                        oldval = datatable_get(otable, simlevel + 2, tmp1, &nullflag);
                        if (nullflag)
                            oldval = 0.0;
                        datatable_set(otable, simlevel + 2, tmp1, oldval + 1.0);
                    }
                    #pragma omp atomic
                    (*p_n)++;
                    if ((*p_n)>nextprint) {
                        n = (*p_n);
                        if(nfiles>100)
                            nextprint += nfiles/100;
                        else
                            nextprint=n;
                        args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files compared\r", n, nfsquare, (int)(((double)n * 100) / (double)nfiles));
                    }

                }
            }
        }
        else {
            for (ifl_c = ifl_added + 1; ifl_c < nfiles;ifl_c++) {
                simlevel = signature_taxsimlevel(signatures[ifl_added], signatures[ifl_c]);
                if (onlylevel < -1 || simlevel == onlylevel) {
                    if (signatures[ifl_added] && signatures[ifl_c])
                        dist = mf((double*)(signatures[ifl_added]->data), (double*)(signatures[ifl_c]->data), sig_len, metharg);
                    else
                        dist = 1.0;
                    if (dist >= 0.0) {
                        tmp1 = (size_t)floor(dist / binwidth);
                        nullflag = 0;
                        oldval = datatable_get(otable, simlevel + 2, tmp1, &nullflag);
                        if (nullflag)
                            oldval = 0.0;
                        datatable_set(otable, simlevel + 2, tmp1, oldval + 1.0);
                    }
                    #pragma omp atomic
                    (*p_n)++;
                    if ((*p_n)>nextprint) {
                        n = (*p_n);
                        if(nfsquare>100)
                            nextprint += nfsquare/100;
                        else
                            nextprint=n;
                        args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) comparisons performed\r", n, nfsquare, (int)(((double)n * 100) / (double)nfsquare));
                    }

                }
            }
        }
    }
    return otable;
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
    int sort_by_distance;
    double* ocmp_results;
    int* ocmp_types;
    datatable_t* otable;
    datatable_t** otable_sub;
    PF_t* outputf;
    size_t i,n;
    size_t* p_n;
    size_t nfsquare;
    size_t countsperproc;
    size_t orderedamount;
    /* parse relevant arguments */
    nocalc = set_method_function(args, &mf, &metharg);
    args_report_info(NULL, "Method function selected successfully\n", numthreads);
    
    outputmode = OUTPUTMODE_NORMAL;
    ocmp_results = NULL;
    ocmp_types = NULL;
    sort_by_distance = args_ispresent(args, "sort_by_distance");
    if (args_ispresent(args, "ordered") || sort_by_distance) {
        outputmode = OUTPUTMODE_ORDERED;
    }
    if (args_ispresent(args, "distancematrix"))
        outputmode = OUTPUTMODE_MATRIX;
    if (args_ispresent(args, "histogram")) {
        outputmode = OUTPUTMODE_HIST;
        binwidth = args_getdouble(args, "histogram", 0, 0.001);
    }
    
    if (outputmode == OUTPUTMODE_ORDERED) {
        args_report_info(NULL, "Output mode: ORDERED\n");
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
        orderedamount = nfiles;
    }
    else {
        issearch = 0;
        orderedamount = nfiles*(nfiles - 1) / 2;
    }
    if (outputmode == OUTPUTMODE_ORDERED) {
        ocmp_results = (double*)malloc(sizeof(double)*orderedamount);
        ocmp_types = (int*)malloc(sizeof(int)*orderedamount);
    }

    PFopen(&outputf, args_getstr(args, "output", 0, "stdout"),"wb");
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
        nfsquare = nfiles*(nfiles - 1) / 2;
        countsperproc = (nfiles / numthreads);
        if (countsperproc*numthreads < nfiles)
            countsperproc++;
        args_report_info(NULL, "Each thread handles " _LLD_ " files\n", countsperproc);
        if (outputmode == OUTPUTMODE_NORMAL) {
            PFprintf(outputf, "File1,File2,Expected_Relation,Distance\n");
        }
        if (outputmode == OUTPUTMODE_HIST) {
            otable_sub = calloc(numthreads, sizeof(datatable_t*));
        }
        else {
            otable_sub = NULL;
        }
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
            size_t ordered_offset=0;
            if (outputmode == OUTPUTMODE_HIST) {
                otable_sub[cthread] = GenDisCal_perform_comparisons_hist_thread(args, signatures, nfiles, sig_len, cthread, numthreads, p_n, mf, metharg);
            }
            else{
            if (firstf > nfiles) firstf = nfiles;
            for (ifl_added = firstf;ifl_added < lastf;ifl_added += numthreads) {
                if (cthread == 0) {
                    args_report_info(NULL, "\t\t\t(" _LLD_ ")", ifl_added);
                }
                if (outputmode == OUTPUTMODE_MATRIX) {
                    datatable_set(otable, ifl_added, ifl_added, 0.0);
                }
                if (issearch) {
                    if (outputmode == OUTPUTMODE_ORDERED) {
                        ocmp_types[ifl_added] = -3;
                    }
                    simlevel = signature_taxsimlevel(signatures[ifl_added], signatures[nfiles]);
                    if (onlylevel < -1 || simlevel == onlylevel) {
                        if (signatures[ifl_added] && signatures[nfiles])
                            dist = mf((double*)(signatures[ifl_added]->data), (double*)(signatures[nfiles]->data), sig_len, metharg);
                        else
                            dist = 1.0;
                        if (dist >= minval && dist <= maxval) {
                            if (outputmode == OUTPUTMODE_NORMAL) {
                                if (signatures[ifl_added] && signatures[ifl_added]->fname) fn1 = signatures[ifl_added]->fname;
                                else fn1 = "<File could not opened>";
                                if (signatures[nfiles] && signatures[nfiles]->fname) fn2 = signatures[nfiles]->fname;
                                else fn2 = "<File could not opened>";
#pragma omp critical
                                PFprintf(outputf, "%s,%s,%s,%f\n", fn1, fn2, taxsimtostr(simlevel), dist);
                            }
                            if (outputmode == OUTPUTMODE_ORDERED) {
                                ocmp_results[ifl_added] = dist;
                                ocmp_types[ifl_added] = simlevel;
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
                            if (nfiles < 100 || (*p_n) % (nfiles / 100) == 0) {
#pragma omp critical
                                args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files compared\r", (*p_n), nfiles, (int)(((double)n * 100) / (double)nfiles));
                            }
                        }
                    }
                }
                else {
                    ordered_offset = 1 + ifl_added*nfiles - (ifl_added + 1)*(ifl_added + 2) / 2;
                    /* note: ordered_offset cannot be negative, but the above expression would return -1 for ifl_added == 0
                     *       this is why "1+" is added at the beginning of the expression and subtracted later.
                     *       While not doing either step should still work in theory as all the variables involved have
                     *       the same size, this approach was preferred to avoid potential problems if this were to change,
                     *       since (uint32_t)(-1) != (uint64_t)(-1), for example.
                     */
                    for (ifl_c = ifl_added + 1; ifl_c < nfiles;ifl_c++) {
                        if (outputmode == OUTPUTMODE_ORDERED)
                            ocmp_types[ifl_c + ordered_offset - 1] = -3;
                        simlevel = signature_taxsimlevel(signatures[ifl_added], signatures[ifl_c]);
                        if (onlylevel < -1 || simlevel == onlylevel) {
                            if (signatures[ifl_added] && signatures[ifl_c])
                                dist = mf((double*)(signatures[ifl_added]->data), (double*)(signatures[ifl_c]->data), sig_len, metharg);
                            else
                                dist = 1.0;
                            if (dist >= minval && dist <= maxval) {
                                if (outputmode == OUTPUTMODE_NORMAL) {
                                    if (signatures[ifl_added] && signatures[ifl_added]->fname) fn1 = signatures[ifl_added]->fname;
                                    else fn1 = "<File could not opened>";
                                    if (signatures[ifl_c] && signatures[ifl_c]->fname) fn2 = signatures[ifl_c]->fname;
                                    else fn2 = "<File could not opened>";
#pragma omp critical
                                    PFprintf(outputf, "%s,%s,%s,%f\n", fn1, fn2, taxsimtostr(simlevel), dist);
                                }
                                if (outputmode == OUTPUTMODE_ORDERED) {
                                    /* see above note for an explanation as to -1 is present here */
                                    ocmp_results[ifl_c + ordered_offset - 1] = dist;
                                    ocmp_types[ifl_c + ordered_offset - 1] = simlevel;
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
            }
        }
        if (outputmode == OUTPUTMODE_HIST) {
            datatable_addtables(otable, otable_sub, numthreads);
            for (i = 0;i < (size_t)numthreads;i++) {
                if (otable_sub[i])datatable_free(otable_sub[i]);
                otable_sub[i] = NULL;
            }
            free(otable_sub);
        }
        args_report_progress(NULL, "A Total of " _LLD_ " comparisons were performed                        \n", n);
    }
    if (!nocalc && otable) {
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
    if (!nocalc && outputmode == OUTPUTMODE_ORDERED) {
        GenDisCal_print_ordered_output(orderedamount, ocmp_types, ocmp_results, nfiles, issearch, sort_by_distance, outputf, signatures);
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
    size_t numfiles, nfilesmod;
    size_t siglen;
    int numthreads;
    filelist = GenDisCal_get_file_list(args, &numfiles);
    if (numfiles < 1) {
        args_report_error(NULL, "Number of files (" _LLD_ ") to be loaded is below 2! Aborting.\n", numfiles);
    }
    else {
        numthreads = args_getint(args, "nprocs", 0, autothreads());
        args_report_info(NULL, "Using up to %d threads\n", numthreads);
        nfilesmod = numfiles;
        siglist = GenDisCal_get_signatures(args, filelist, &nfilesmod, &siglen, numthreads);
        if (nfilesmod < 1) {
            args_report_error(NULL, "Number of files (" _LLD_ ") to be loaded is below 2! Aborting.\n", numfiles);
        }
        else {
            if(!args_ispresent(args,"keepsigs") || strcmp(args_getstr(args,"keepsigs",0,"continue"),"stop")!=0)
                GenDisCal_perform_comparisons(args, siglist, nfilesmod, siglen, numthreads);

            for (i = 0;i < numfiles;i++) {
                if (siglist[i])signature_free(siglist[i]);
                siglist[i] = NULL;
            }
        }
        if (siglist)
            free(siglist);
    }
    for (i = 0;i < numfiles;i++) {
        free(filelist[i]);
    }
    free(filelist);
    return 0;
}

#if TEST
#else
int main(int argc, char** argv) {
    args_t* args;
    int result;
    args = GenDisCal_init_args(argc, argv);
    if(!args_is_helpmode(args))
        result = GenDisCal(args);
    args_free(args);
#ifdef _DEBUG
    system("PAUSE");
#endif
    return result;
}
#endif