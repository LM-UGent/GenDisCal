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

/* GLOBAL DEFINES */
#define VERSION_NAME    "GenDisCal v1.2.2"
#define CHANGES \
"v1.2.2\n"\
"Should no longer crash in cases where signatures cannot be computed.\n"\
"A warning will now be issued when supplying protein fasta instead instead of nucleotide fasta.\n"\
"v1.2.1\n"\
"Corrected analysis of blast output.\n"\
"Lowercase nucleotide fasta should now be handled properly.\n"\
"v1.2\n" \
"-s / --search now works with .sdb files" \
"v1.1\n" \
"The program has been fully reworked internally, with hopefully not consequences on normal functionality.\n" \
"In addition, the following options have been added:\n" \
"--license         displays license\n" \
"--version         displays current version\n" \
"--changes         displays changes since the last minor update\n" \
"-c / --canonical  shorten signatures to only include canonical k-mers\n" \
"-a / --filealias  add aliasias to be displayed in place of filenames in the output \n" \
"-e / --exclude    exclude samples based on taxonomy \n" \
"-i / --include    include samples based on taxonomy \n" \
"-u / --clustering perform clustering based on the distance matrix"
#ifdef _DEBUG
/*#define NOOPENMP*/
#endif
/* INCLUDES */
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
#include "genosig.h"
#include "textparsing.h"
#include "taxonomy.h"
#include "datatable.h"
#include "suffixtree.h"

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
    args_add(result, "license", 0, "");
    args_add(result, "version", 0, "");
    args_add(result, "changes", 0, "");
    args_add(result, "basis", 'b', "str,var");
    args_add(result, "method", 'm', "str,var");
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
    args_add(result, "canoncial", 'c', "");
    args_add(result, "filealias", 'a', "str");
    args_add(result, "exclude", 'e', "...");
    args_add(result, "include", 'i', "...");
    args_add(result, "clustering", 'u', "str");
    args_add(result, "z", 'z', "");
#ifndef NOOPENMP
    args_add(result, "nprocs", 'n', "int");
#endif
    args_add_help(result, "license", "LICENSE", "Prints the license(s) associated with GenDisCal.\n", "Print license(s)");
    args_add_help(result, "version", "LICENSE", "Prints the license(s) associated with GenDisCal.\n", "Print license(s)");
    args_add_help(result, "changes", "LICENSE", "Prints the license(s) associated with GenDisCal.\n", "Print license(s)");
    args_add_help(result, "basis", "SIGNATURE BASIS",
        "Signature bases are used to transform a genome into a vector of numbers. "
        "Available signatures are as follows:\n"
        "   len         -- use genome length as a signature\n"
        "   GC          -- use GC content as a signature\n"
        "   freq <k>    -- use frequencies of <k>-mers as a signature\n"
        "   karl <k>    -- use karlin signatures of <k>-mers\n"
        "   exp <k>     -- use <k>-mer bias with regards to monomer composition\n"
        /*"   gene <name> -- use nucleotide identity for the gene stored in the file <name>\n"
        "                  in addition, <name> may be \"16S\".\n"
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
        "   reldist        -- use the mean of pairwise relative distances between signature values\n"
        "   ANIb <path>    -- use nucleotide identity computed by <path> (full path to blastn)\n",
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
        " > approxANI (no equivalent)\n"
        "    approxANI is a minhash-based approximation of Average Nucleotide\n"
        "    Identity (ANI). This method is similar to the one used by the \n"
        "    Mash software (Ondov et al., Gen. biol., 2016)\n"
        "    The ambiguity region for species is [0.04-0.06]\n"
        " > combinedSpecies (no equivalent)\n"
        "    combinedSpecies produces a probability that two genomes belong to\n"
        "    the same bacterial species, based on multiple metrics\n"
        /*" > 16S (-b gene 16S -m alignseqs)"
        "    Nucleotide identity of aligned 16S rRNA genes is usually considered to\n"
        "    be a good indication of two bacterial genomes belonging to the same\n"
        "    species.\n"
        "    The typically used thresholds are 0.0135 and 0.03, with values below\n"
        "    the threshold belonging to the same species. We note that in our tests,\n"
        "    2.5%% of <same genus> comparisons fell below 0.0135, and 18.6%% fell below\n"
        "    0.03.\n"*/,
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
        "  GTDB    -- filename[tab]d__Domain;p__Phylum;c__Class;...\n"
        "  NCBI    -- filename,NCBI_taxid\n"
        "When using the NCBI format, please ensure that nodes.dmp.txt is present in the folder where "
        "the command is run.\n"
        "If the the supplied filename ends with '.ttd', it is assumed to be a pre-computed taxonomy "
        "file, and the second argument does not need to be specified. Otherwise, such a file will be "
        "generated by this program.\n",
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
        "The output format is a always a comma-separated file.\n"
        "This option overrides --distancematrix and --clustering\n",
        "enable histogram output");
    args_add_help(result, "distancematrix", "DISTANCE MATRIX OUTPUT",
        "If this option is present, the output will be a distance matrix instead of a list of "
        "comparisons.\n"
        "The output format is a always a comma-separated file.\n"
        "This option is required for --clustering\n",
        "enable distance matrix output");
    args_add_help(result, "clustering", "CLUSTERING OUTPUT",
        "If this option is present, a tree will be produced based on the specified method.\n"
        "The output format is a newick tree.\n"
        "available methods: ward (default)\n",
        "enable clustering output");
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
        "  Same_Subgroup\n"
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
    args_add_help(result, "canoncial", "USE CANONICAL SIGNATURES",
        "If this option is present, elements of k-mer-based signatures which to the reverse "
        "complement of 'smaller' k-mers signatures will be removed. This can accelerate calculations.\n"
        "We consider that k-mer can be compared in the same way as numbers where each letter "
        "represents a digit, with A<C<G<T. For example, under this notation AAAG < CTTT, and hence "
        "only the value associated with AAAG will be kept.\n",
        "Simplify signatures to keep canoncial elements only");
    args_add_help(result, "filealias", "USE FILE ALIASES IN OUTPUT",
        "If this option is present, elements in the output will be renamed according to taxonomy "
        "or aliases contained in the alias file supplied with this option.\n",
        "Change the name of sequences in the output");
    args_add_help(result, "exclude", "EXCLUDE SPECIFIC TAXONOMIC NAMES",
        "All names following this option will not be included for comparison. Any genome whose "
        "taxonomy contains any of these names will therefore be excluded. This remains true "
        "even if these sequences match one of the names specified with -i/--include.\n",
        "Exclude taxids");
    args_add_help(result, "include", "INCLUDE SPECIFIC TAXONOMIC NAMES",
        "When this flag is present, only the sequences which contains at least one of the names "
        "specified after this option will be included. Sequences included this way may still "
        "be excluded using -e/--exclude.\n",
        "Include taxids");
    args_add_help(result, "z", "END OF OPTION",
        "Use this flag to terminate arguments with unspecified counts. All names specified after "
        "this flag (but before the next flag/option) will be considered as query file names.\n",
        "end of option");
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
            PFclose(f);
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

int set_basis_function(args_t* args, genosig_bf* bf, size_t* sig_len, size_t* kmer_len) {
    size_t maxi, i;
    size_t kmerlen;
    int nocalc;
    char* basisstr;
    char* presetstr;
    char* possible_bases[] = { "GC","freq","karl","karlL","exp","mmzinf","len","minhash","sequence","findgene","filename","none" };
    genosig_bf bf_list[] = { genosig_GC,genosig_kmerfreq,genosig_karlinsig,genosig_karlinsigL,genosig_mmz,genosig_mmz_inf, genosig_length,genosig_minhash,NULL,NULL,genosig_keepfilenameonly,genosig_keepfilenameonly };
    nocalc = 0;
    presetstr = args_getstr(args, "preset", 0, NULL);
    if (presetstr) {
        if (strcmp(presetstr, "GC") == 0) {
            *bf = genosig_GC;
            kmerlen = 1;
            *sig_len = 1;
        }
        if (strcmp(presetstr, "TETRA") == 0) {
            *bf = genosig_mmz;
            kmerlen = 4;
            *sig_len = (1LL << (4 * 2));
        }
        else if (strcmp(presetstr, "PaSiT4") == 0) {
            *bf = genosig_karlinsig;
            kmerlen = 4;
            *sig_len = (1LL << (4 * 2));
        }
        else if (strcmp(presetstr, "PaSiT6") == 0) {
            *bf = genosig_karlinsig;
            kmerlen = 6;
            *sig_len = (1LL << (6 * 2));
        }
        else if (strcmp(presetstr, "approxANI") == 0) {
            *bf = genosig_minhash;
            kmerlen = 21;
            *sig_len = 3000;
        }else {
            args_report_warning(NULL, "Unknown preset <%s>, using default basis: <karl 4>\n", presetstr);
            kmerlen = 4;
            *bf = genosig_karlinsig;
        }
    }
    else {
        basisstr = args_getstr(args, "basis", 0, "karl");
        if (basisstr[0] == 's' && basisstr[1] == '_') basisstr += 2;
        maxi = sizeof(possible_bases) / sizeof(char*);
        for (i = 0;i < maxi;i++) {
            if (strcmp(possible_bases[i], basisstr) == 0)break;
        }
        if (i >= maxi) {
            if (beginswith("karl", basisstr)) {
                if (basisstr[4] == '*') {
                    kmerlen = atoll(basisstr + 5);
                    *bf = genosig_karlinsig;
                    args_report_warning(NULL, "'-b karl*%d' is deprecated, use '-b karl %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
                else {
                    kmerlen = atoll(basisstr + 4);
                    *bf = genosig_karlinsig;
                    args_report_warning(NULL, "'-b karl%d' is deprecated, use '-b s_karl %d' instead\n",(int)kmerlen, (int)kmerlen);
                }
            }
            else if (beginswith("freq", basisstr)) {
                if(basisstr[4] == '*') {
                    kmerlen = atoll(basisstr + 5);
                    *bf = genosig_kmerfreq;
                    args_report_warning(NULL, "'-b freq*%d' is deprecated, use '-b freq %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
                else {
                    kmerlen = atoll(basisstr + 4);
                    *bf = genosig_kmerfreq;
                    args_report_warning(NULL, "'-b freq%d' is deprecated, use '-b s_freq %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
            }
            else if (beginswith("markz", basisstr)) {
                if (basisstr[5] == '*') {
                    kmerlen = atoll(basisstr + 6);
                    *bf = genosig_mmz;
                    args_report_warning(NULL, "'-b markz*%d' is deprecated, use '-b exp %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
                else {
                    kmerlen = atoll(basisstr + 5);
                    *bf = genosig_mmz;
                    args_report_warning(NULL, "'-b markz%d' is deprecated, use '-b s_exp %d' instead\n", (int)kmerlen, (int)kmerlen);
                }
            }
            else if (strcmp("TETRA", basisstr) == 0) {
                *bf = genosig_mmz;
                kmerlen = 4;
                args_report_warning(NULL, "'-b TETRA' is deprecated, use '-p TETRA' instead\n");
            }
            else {
                args_report_warning(NULL, "Unknown basis <%s>, defaulting to <karl>\n", basisstr);
                *bf = genosig_karlinsig;
                basisstr = "karl";
            }
        }
        else {
            *bf = bf_list[i];
            if (*bf == genosig_minhash) {
                kmerlen = 21;
                *sig_len = 3000;
            }
            else {
                kmerlen = 4;
            }
        }
        if (strcmp(basisstr, "GC") == 0) {
            args_report_info(NULL, "Using GC content as signatures\n");
        }
        else if (strcmp(basisstr, "len") == 0) {
            args_report_info(NULL, "Using genome lengths as signatures\n");
        }
        else {
            kmerlen = atoll(args_getstr(args, "basis", 1, "4"));
            if (kmerlen == 0) kmerlen = 4;
            args_report_info(NULL, "Calculating %s signatures\n", basisstr);
            args_report_info(NULL, "k-mer length set to %d (may not be meaningful)\n", (int)kmerlen);
        }
        if (kmerlen < 0) kmerlen = 0;
    }
    *kmer_len = kmerlen;
    return nocalc;
}
int set_method_function(args_t* args, genosig_df* mf, any_t* metharg) {
    char* methstr;
    char* presetstr;
    any_t resarg;
    presetstr = args_getstr(args, "preset", 0, NULL);
    resarg.d = 0.02;
    if (presetstr) {
        if (strcmp(presetstr, "TETRA") == 0) {
            *mf = genodist_pearscorr;
            resarg.d = 0.0;
        }
        else if (strcmp(presetstr, "PaSiT4") == 0) {
            *mf = genodist_satman;
            resarg.d = 0.02;
        }
        else if (strcmp(presetstr, "PaSiT6") == 0) {
            *mf = genodist_satman;
            resarg.d = 0.02;
        }
        else if (strcmp(presetstr, "approxANI") == 0) {
            *mf = genodist_approxANI;
            resarg.d = 21.0;
        }
        else {
            args_report_warning(NULL, "Unknown preset <%s>, using default method: <PaSiT 0.02>\n", presetstr);
            *mf = genodist_satman;
            resarg.d = 0.02;
        }
    }
    else {
        methstr = args_getstr(args, "method", 0, "SVC");
        if (strcmp(methstr, "AMD") == 0 || strcmp(methstr, "manhattan") == 0) {
            *mf = genodist_manhattan;
            resarg.d = 0.0;
        }
        else if (strcmp(methstr, "ED") == 0 || strcmp(methstr, "euclid") == 0) {
            *mf = genodist_euclidian;
            resarg.d = 0.0;
        }
        else if (strcmp(methstr, "SVC") == 0 || strcmp(methstr, "PaSiT") == 0 || strcmp(methstr, "satman") == 0) {
            *mf = genodist_satman;
            resarg.d = args_getdoublefromstr(args, "method", 1, 0.02);
        }
        else if (strcmp(methstr, "PaSiTL") == 0) {
            *mf = genodist_PaSiTL;
            resarg.d = args_getdoublefromstr(args, "method", 1, 0.02);
        }
        else if (strcmp(methstr, "sateucl") == 0 ) {
            *mf = genodist_sateucl;
            resarg.d = args_getdoublefromstr(args, "method", 1, 0.02);
        }
        else if (strcmp(methstr, "hamming") == 0) {
            *mf = genodist_satman;
            resarg.d = 0.0;
        }
        else if (strcmp(methstr, "corr") == 0 || strcmp(methstr, "pearson") == 0) {
            *mf = genodist_pearscorr;
            resarg.d = 0.0;
        }
        else if (strcmp(methstr, "rankcorr") == 0 || strcmp(methstr, "spearman") == 0) {
            *mf = genodist_rankcorr;
            resarg.d = 0.0;
        }
        else if (strcmp(methstr, "reldist") == 0) {
            *mf = genodist_manhattan;
            resarg.d = 0.0;
        }
        else if (beginswith("SVC", methstr)) {
            resarg.d = atof(methstr + 3);
            *mf = genodist_satman;
            args_report_warning(NULL, "'-m SVC<f>' is deprecated, use '-m SVC <f>' instead\n");
        }
        else if (beginswith("alignseqs", methstr)) {
            resarg.d = 1;
            *mf = genodist_ANI;
        }
        else if (beginswith("ANIb", methstr)) {
            resarg.str = args_getstr(args,"method",1,NULL);
            if (resarg.str == NULL) args_report_error(NULL, "Please specify the path to blast to use ANIb\n");
            *mf = genodist_externANIb;
        }
        else if (strcmp(methstr, "approxANI") == 0) {
            *mf = genodist_approxANI;
            resarg.d = args_getdoublefromstr(args, "method", 1, 0.02);
            resarg.d = 21.0;
        }
        else {
            resarg.str = args_getstr(args, "method", 1, NULL);
            if (resarg.str != NULL) {
                args_report_warning(NULL, "Unknown method <%s %s>, using default: <PaSiT 0.02>\n", presetstr);
            }
            else {
                args_report_warning(NULL, "Unknown method <%s>, using default: <PaSiT 0.02>\n", presetstr);
            }
            *mf = genodist_satman;
        }
    }
    *metharg = resarg;
    return 0;
}

void write_sigfile(const char* basename, genosig_t* source) {
    PF_t* f;
    char* newname;
    size_t baselen;
    baselen = strlen(basename);
    newname = (char*)malloc(baselen + 5);
    memcpy(newname, basename, baselen);
    memcpy(newname + baselen, ".sig", 5);
    PFopen(&f, newname, "wb");
    genosig_savebin(source,f);
    PFclose(f);
}
genosig_t* get_signatures(char* filename, genosig_bf basis, size_t nlen, size_t minseqlen, int nonredundant, int singlestrand, int keepgenome, genosig_t* refgenes, size_t winsize) {
    genosig_t* result;
    nucseq** nsarr;
    size_t outcount;
    size_t badcount;
    PF_t* f;
    if (endswith(".sig", filename) || endswith(".sdb", filename)) {
        result = genosig_import(filename);
    }
    else {
        result = NULL;
        if (basis != genosig_keepfilenameonly) {
            PFopen(&f, filename, "rb");
            if (f) {
                nsarr = nucseq_array_from_fasta(f, &outcount, 1, minseqlen,&badcount);
                if (badcount) {
                    args_report_warning(NULL, "File %s contained sequences which may not be nucleotide fasta.\n", filename);
                }
                if (outcount == 0) {
                    args_report_warning(NULL, "File %s did not contain any sequences longer than " _LLD_ " nucleotides.\n", filename, minseqlen);
                }
                result = genosig_fullgenome(nsarr, outcount, singlestrand, COPYLVL_INTEGRATE);
                if (winsize > 0) result = genosig_makewindows(result, winsize, 1);
                if (!genosig_info_name(result))
                    result = genosig_autoname(result, filename);
                if (refgenes) {
                    result = genosig_findgenes(result, singlestrand, refgenes, 0);
                }
                if (basis != NULL) {
                    result = basis(result, nlen);
                    if (nonredundant) result = genosig_keepnonredundant(result, nlen);
                }
                else if (keepgenome) {
                    result = genosig_suffixtrees(result);
                }
                if (!keepgenome) {
                    genosig_cleargenomedata(result);
                }
                PFclose(f);
            }
        }
        else {
            result = genosig_genomefileref(filename);
            result = genosig_keepfilenameonly(result,0);
        }
    }
    return result;
}

genosig_t* get_gene_signature(char* gene, genosig_bf basis, size_t nlen, int nonredundant, int singlestrand, int keepgenome, size_t winsize) {
    genosig_t* result;
    nucseq** nsarr;

    char* SSU = 
        "TAATTGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGTATTGAAAGGTGCTTGCACCTGGACGAGTGGCGGACGGGTGAGTAACACGTGGGAACCTGCCCTT"
        "AAGTGGGGGATAACATTTGGAAACAGATGCTAATACCGCATAAAACCGCATGGTTAAAGGCTGAAAGTGGCGTGAGCTATCGCTTTTGGATGGGCCCGCGTCGGATTAGCTAGTTGGTGAGGTAACGGCTCAC"
        "CAAGGCGACGATCCGTAGCCGGTCTGAGAGGATGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCAATGCCGCG"
        "TGAGTGAAGAAGGCCTTCGGGTTGTAAAGCTCTTTTGTTGGGGAAGAAAGGTCGGCAGGTAACTGTTGTCGGCGTGACGGTACCCAACGAAAAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGT"
        "AGGGTGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGAGCGCAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCTCGGCTTAACCGGGGAAGTGCATTGGAAACTGGGAAACTTGAGTGCAGAAGAGGAGAG"
        "TGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGT"
        "CCACGCCGTAAACGATGAATGCTAGGTGTTGGGGGGTTTCCGCCCTTCAGTGCCGCAGCTAACGCATTAAGCATTCCGCCTGGGGAGTACGGCCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCAC"
        "AAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCTTTTGAACACCTTAGAGATAAGGTTTTCCCTTCGGGGACAAAATGACAGGTGCTGCATGGTTGTCGTCAGC"
        "TCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCATTAGTTGCCAGCATTAAGTTGGGCACTCTAGTGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATC"
        "ATGGCCCTTATGACCTGGGCTACACACGTGCTACAATGGATGGTACAAAGAGTTGCGAGACCGCGAGGTCAAGCTAATCTCTTAAAGCCATTCTCAGTTCGGATTGTAGTCTGCAACTCGACTACATGAAGTC"
        "GGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTTGTAACACCCGAAGTCGGTGGCCTAACCTTAGGGAGGGAGCCGACTAA"
        "GGTGGGACAGATGACTGGGGTGAAGTCGTAACAAGGTAGCCGTAGGGGAACCTGCGGCTGGATCACCTCCTTTCT";

    nsarr = (nucseq**)malloc(sizeof(nucseq*));
    nsarr[0] = (nucseq*) calloc(1, sizeof(nucseq));
    if(strcmp(gene,"SSU")==0 || strcmp(gene,"16S")==0)
        nucseq_from_string(nsarr[0], SSU);
    result = genosig_fullgenome(nsarr, 1, singlestrand, COPYLVL_INTEGRATE);
    if (winsize > 0) result = genosig_makewindows(result, winsize, 1);
    if (basis != NULL) {
        result = basis(result, nlen);
        if (nonredundant) result = genosig_keepnonredundant(result, nlen);
    }
    else if (keepgenome) {
        result = genosig_suffixtrees(result);
    }
    if (!keepgenome) {
        genosig_cleargenomedata(result);
    }
    return result;
}


genosig_t** GenDisCal_get_signatures(args_t* args, char** sourcefiles, size_t* numfiles, size_t* sig_len, int numthreads, size_t* p_nsearch) {
    /* This function is too massive and should probably be split into its consitituent parts */
    genosig_bf bf;
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
    taxtree_t* taxtree;
    char* line;
    char* taxtype;
    char* badload;
    char* badtax;
    char* tmpstr;
    char* tmpstr2;
    PF_t* f;
    PF_t* fextra;
    size_t taxfilenamelen;
    genosig_t** result;
    genosig_t** tmpresult;
    genosig_t** subres;
    genosig_t* refgenes;
    genosig_t* allsigs;
    genosig_t* sigavg;
    int gensigs;
    int nonred;
    int singlestrand;
    int usealiases;
    double renaminglevelmin;
    double renaminglevelmax;
    size_t nfiles;
    size_t winsize;
    char** include_list;
    char** exclude_list;
    int64_t* include_ids;
    int64_t* exclude_ids;
    size_t nincl;
    size_t nexcl;
    time_t step_start;
    time_t step_duration;
    int hours, minutes, seconds;
    int shouldbeincluded;
    
    renaminglevelmin = 6.0;
    renaminglevelmax = 7.0;
    nfiles = *numfiles;
    nocalc = set_basis_function(args, &bf, sig_len, &kmerlen);
    winsize = 0;
    if (bf == genosig_minhash) {
        winsize = kmerlen;
        kmerlen = *sig_len;
    }

    result = (genosig_t**)calloc(nfiles, sizeof(genosig_t*));
    if (!nocalc) {
        /* prepare taxonomy */
        taxtree = NULL;
        taxtype = NULL;
        usealiases = 0;
        if (args_ispresent(args, "taxonomy")) {
            step_start = time(NULL);
            line = args_getstr(args, "taxonomy", 0, "taxonomy.csv");
            taxtype = args_getstr(args, "taxonomy", 1, "PCOF_S");
            if (PFopen(&f, line, "rb")) args_report_warning(NULL, "Taxonomy file <%s> could not be loaded\n", line);
            else {
                args_report_progress(NULL, "Loading taxonomy...\n");
                if (endswith(".ttd", line)) {
                    taxtree = taxtree_load(f);
                }
                else {
                    if (strcmp(taxtype, "PCOF_S") == 0) {
                        taxtree = taxtree_fromlineages(f, TAXFMT_PCOF_S);
                    }
                    else if (strcmp(taxtype, "GTDB") == 0) {
                        taxtree = taxtree_fromlineages(f, TAXFMT_GTDB);
                    }
                    else if (strcmp(taxtype, "NCBI") == 0) {
                        tmpstr = os_extractdirname(line);
                        n = strlen(tmpstr);
                        k = strlen("nodes.dmp.txt");
                        tmpstr2 = (char*)malloc(n + k + 1);
                        memcpy(tmpstr2, tmpstr, n);
                        memcpy(tmpstr2+n, "nodes.dmp.txt", k + 1);
                        if (PFopen(&fextra, tmpstr2, "rb")) {
                            args_report_warning(NULL, "<nodes.dmp.txt> could not be found. Taxonomy not loaded.\n");
                            args_report_info(NULL, "Please add the file <nodes.dmp.txt> obtained from the NCBI FTP to the current directory\n");
                        }
                        else {
                            taxtree = taxtree_fromrelations(fextra, TAXFMT_NCBI);
                            PFclose(fextra);
                            k = strlen("merged.dmp.txt");
                            tmpstr2 = (char*)realloc(tmpstr2, n + k + 1);
                            memcpy(tmpstr2 + n, "merged.dmp.txt", k + 1);
                            if (PFopen(&fextra, tmpstr2, "rb")) {
                                args_report_info(NULL, "<merged.dmp.txt> not found. Some taxids may be incorrect.\n", line);
                            }
                            else {
                                taxtree_rename(taxtree, fextra, TAXFMT_NCBI);
                                PFclose(fextra);
                            }
                            k = strlen("names.dmp.txt");
                            tmpstr2 = (char*)realloc(tmpstr2, n + k + 1);
                            memcpy(tmpstr2 + n, "names.dmp.txt", k + 1);
                            if (PFopen(&fextra, tmpstr2, "rb")) {
                                args_report_info(NULL, "<nmes.dmp.txt> not found. Readable names will not be understood.\n", line);
                            }
                            else {
                                taxtree_addleaves(taxtree, fextra, TAXFMT_NCBI);
                                PFclose(fextra);
                            }
                            taxtree_addleaves(taxtree, f, TAXFMT_CSV);
                        }
                        free(tmpstr);
                        free(tmpstr2);
                    }
                    if (taxtree) {
                        taxfilenamelen = strlen(line);
                        tmpstr = malloc(taxfilenamelen + 5);
                        memcpy(tmpstr, line, taxfilenamelen);
                        memcpy(tmpstr + taxfilenamelen, ".ttd", 5);
                        if (PFopen(&fextra, tmpstr, "wb")) {
                            args_report_warning(NULL, "<%s> could not be created. Please make sure you have write permission in this folder.\n", tmpstr);
                        }
                        else {
                            taxtree_save(taxtree, fextra);
                            PFclose(fextra);
                        }
                        free(tmpstr);
                    }
                }
                PFclose(f);
                step_duration = time(NULL) - step_start;
                hours = (int)(step_duration / (60 * 60));
                minutes = (int)((step_duration % (60 * 60)) / (60));
                seconds = (int)((step_duration % (60)));
                args_report_progress(NULL, "Taxonomy loaded (%d hr %d min %d sec)\n", hours, minutes, seconds);
            }
        }
        if (args_ispresent(args, "filealias")) {
            tmpstr = args_getstr(args, "filealias", 0, NULL);
            if (taxtree || tmpstr) {
                if (!taxtree) taxtree = taxtree_alloc();
                else usealiases = RENAMEBY_SPECIES;
                if (tmpstr) {
                    if (PFopen(&f, tmpstr, "r")) {
                        args_report_warning(NULL, "Alias-containing file <%s> could not be opened\n", tmpstr);
                    }
                    else{
                        taxtree_rename(taxtree, f, TAXFMT_CSV);
                        PFclose(f);
                        usealiases = RENAMEBY_NODELABEL;
                    }
                }
            }
            else {
                args_report_warning(NULL, "Please specify either an alias file using -a <path> or provide taxonomy with -t <path> to use aliases\n");
            }
        }
        /* read the signatures */
        nonred = args_ispresent(args, "nonredundant");
        singlestrand = 0;
        refgenes = NULL;
        if (args_ispresent(args, "basis")) {
            tmpstr = args_getstr(args, "basis", 0, "--");
            if (beginswith("s_", tmpstr))
                singlestrand = 1;
            if (strcmp(tmpstr, "findgene") == 0) {
                tmpstr = args_getstr(args, "basis", 1, "gene.fasta");
                if (strcmp(tmpstr, "16S") == 0 || strcmp(tmpstr, "SSU") == 0) {
                    refgenes = get_gene_signature(tmpstr, NULL, 0, 0, 1, 1, winsize);
                }
                else {
                    refgenes = get_signatures(tmpstr, NULL, 0, 0, 0, 1, 1, NULL,winsize);
                }
                refgenes = genosig_suffixtrees(refgenes);
            }
        }
        gensigs = (args_ispresent(args, "keepsigs")?1:0);
        if (gensigs && args_ispresent(args, "filelist"))gensigs = 2;
        args_report_progress(NULL, "Reading signatures...\n");
        minseqlen = 2000;
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
        args_report_info(NULL, "Signature length is " _LLD_ "\n", kmerlen);
        i = 0;
        i_p = &i;
        step_start = time(NULL);
        #pragma omp parallel
        {
            int nullflag;
            int muxed;
            size_t ifl_added = 0;
            size_t next_ifla_report = 0;
            int cthread = omp_get_thread_num();

            #pragma omp critical
            ifl_added = ((*i_p)++);
            while(ifl_added<nfiles) {
                result[ifl_added] = get_signatures(sourcefiles[ifl_added], bf, kmerlen, minseqlen, nonred, singlestrand, 0, refgenes, winsize);
                muxed = genosig_info_ismuxed(result[ifl_added]);
                if (result[ifl_added]) {
                    if (!muxed && gensigs == 1 && !endswith(".sig", sourcefiles[ifl_added]))
                        write_sigfile(sourcefiles[ifl_added], result[ifl_added]);
                    nullflag = 0;
                    if (taxtree)
                        result[ifl_added] = genosig_importlineage(result[ifl_added], taxtree,&nullflag);
                    if (nullflag) {
                        #pragma omp atomic
                        (*p_notax)++;
                        badtax = sourcefiles[ifl_added];
                    }
                    #pragma omp atomic
                    (*p_n)++;
                }
                else {
                    badload = sourcefiles[ifl_added];
                }
                if (cthread == 0 && (nfiles < 100 || (*p_n) > next_ifla_report)) {
                    args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files loaded\r", n, nfiles, (int)((n * 100) / nfiles));
                    while (next_ifla_report < (*p_n))
                        next_ifla_report += (nfiles / 100) + 1;
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
        /* demultiplex multiplexed signatures*/
        nfiles = n;
        j = 0;
        for (i = 0;i<n;i++) {
            if (result[i]) {
                if (genosig_info_ismuxed(result[i])) {
                    nfiles += genosig_info_length(result[i]) - 1;
                    *p_nsearch = genosig_info_length(result[i]);
                    j++;
                }
                else {
                    *p_nsearch = 1;
                }
            }
        }
        if (!args_ispresent(args, "search")) {
            *p_nsearch = 0;
        }
        if (n != nfiles) {
            args_report_progress(NULL, "Additionally, " _LLD_ " signature(s) are contained inside of " _LLD_ " SDB file(s)\n", nfiles - n + j, j);
            tmpresult = (genosig_t**)malloc(sizeof(genosig_t*)*nfiles);
            extrafiles = 0;
            for (i = 0;i < n;i++) {
                if (result[i] && genosig_info_ismuxed(result[i])) {
                    subres = genosig_demux(result[i],&j);
                    /* rem: demultiplexing invalidates result[i], so it should not be free'd. Insead it should simply be set to NULL.*/
                    result[i] = NULL;
                    memcpy(tmpresult + extrafiles, subres, sizeof(genosig_t*)*(j));
                    free(subres);
                    if (taxtree) {
                        /* assign taxonomy */
                        for (k = 0;k < j;k++) {
                            nflag = 0;
                            if (taxtree) tmpresult[k + extrafiles] = genosig_importlineage(tmpresult[k + extrafiles], taxtree, &nflag);
                            if (usealiases) tmpresult[k + extrafiles] = genosig_renameby(tmpresult[k + extrafiles], taxtree, COPYLVL_FULL, usealiases);
                            if (nflag) {
                                notax++;
                                badtax = sourcefiles[i];
                            }
                        }
                    }
                    extrafiles += j;
                }
                else {
                    tmpresult[extrafiles] = result[i];
                    if (taxtree) tmpresult[extrafiles] = genosig_importlineage(tmpresult[extrafiles], taxtree, &nflag);
                    if (usealiases) tmpresult[extrafiles] = genosig_renameby(tmpresult[extrafiles], taxtree, COPYLVL_FULL, usealiases);
                    extrafiles++;
                }
            }
            free(result);
            result = tmpresult;
            nfiles = extrafiles;
        }
        else {
            for (i = 0;i < n;i++) {
                if (result[i]) {
                    if (taxtree) result[i] = genosig_importlineage(result[i], taxtree, &nflag);
                    if (usealiases) result[i] = genosig_renameby(result[i], taxtree, COPYLVL_FULL, usealiases);
                }
            }
        }
        /* report on taxonomy assignment */
        if (notax > 0) {
            if (notax != 1)
                args_report_warning(NULL, "Taxonomy could not be assigned for " _LLD_ " signatures\n", notax);
            else if (badtax)
                args_report_warning(NULL, "Taxonomy could not be assigned for <%s>\n", badtax);
        }
        step_duration = time(NULL) - step_start;
        hours = (int)(step_duration / (60 * 60));
        minutes = (int)((step_duration % (60 * 60)) / (60));
        seconds = (int)((step_duration % (60)));
        args_report_progress(NULL, "Finished reading signatures (%d hr %d min %d sec)\n", hours, minutes, seconds);
        /* include/exclude signatures */
        include_list = NULL;
        exclude_list = NULL;
        include_ids = NULL;
        exclude_ids = NULL;
        nincl = 0;
        nexcl = 0;
        if (args_ispresent(args, "include")) {
            nincl = args_countargs(args, "include");
            if (args_ispresent(args, "search")) {
                include_list = (char**)malloc(sizeof(char*)*(nincl+1));
                include_ids = (int64_t*)malloc(sizeof(int64_t)*(nincl+1));
                include_list[nincl] = args_getstr(args, "include", (int)nincl, "<none>");
                if (taxtree)
                    include_ids[nincl] = taxtree_getid(taxtree, include_list[nincl]);
                for (i = 0;i < nincl;i++) {
                    include_list[i] = args_getstr(args, "include", (int)i, "<none>");
                    if (taxtree)
                        include_ids[i] = taxtree_getid(taxtree, include_list[i]);
                }
                nincl++;
            }
            else {
                include_list = (char**)malloc(sizeof(char*)*nincl);
                include_ids = (int64_t*)malloc(sizeof(int64_t)*nincl);
                for (i = 0;i < nincl;i++) {
                    include_list[i] = args_getstr(args, "include", (int)i, "<none>");
                    if (taxtree)
                        include_ids[i] = taxtree_getid(taxtree, include_list[i]);
                }
            }
            
        }
        if (args_ispresent(args, "exclude")) {
            nexcl = args_countargs(args, "exclude");
            exclude_list = malloc(sizeof(char*)*nexcl);
            exclude_ids = malloc(sizeof(int64_t)*nexcl);
            for (i = 0;i < nexcl;i++) {
                exclude_list[i] = args_getstr(args, "exclude", (int)i, "<none>");
                if (taxtree)
                    exclude_ids[i] = taxtree_getid(taxtree, exclude_list[i]);
            }
        }
        j = 0;
        for (i = 0;i < nfiles;i++) {
            if (!include_list) {
                shouldbeincluded = 1;
            }
            else {
                shouldbeincluded = 0;
                for (k = 0;k < nincl;k++) {
                    if (include_ids[k] > 0 && genosig_matchesname(result[i], include_ids[k], include_list[k])) {
                        shouldbeincluded = 1;
                        break;
                    }
                }
            }
            if (exclude_list) {
                for (k = 0;k < nexcl;k++) {
                    if (exclude_ids[k] > 0 && genosig_matchesname(result[i], exclude_ids[k], exclude_list[k])) {
                        shouldbeincluded = 0;
                        break;
                    }
                }
            }
            if (shouldbeincluded) {
                result[j] = result[i];
                j++;
            }
            else {
                genosig_free(result[i]);
                result[i] = NULL;
            }
        }
        result = (genosig_t**)realloc(result, sizeof(genosig_t*)*j);
        nfiles = j;
        if (i != j) {
            args_report_progress(NULL, _LLD_ " signatures were kept\n", nfiles);
        }
        if (include_ids) free(include_ids);
        if (include_list) free(include_list);
        if (exclude_ids) free(exclude_ids);
        if (exclude_list) free(exclude_list);
        /* generate signature file */
        if (gensigs) {
            line = args_getstr(args, "filelist", 0, "all_sigs.list");
            if (strcmp(args_getstr(args, "keepsigs", 0, "continue"), "stop") == 0)
                line = args_getstr(args, "output", 0, "all_sigs.list");
            if (strcmp(args_getstr(args, "keepsigs", 0, "continue"), "average") == 0)
                line = args_getstr(args, "output", 0, "all_sigs.list");
            if (strcmp(args_getstr(args, "keepsigs", 0, "continue"), "text") == 0)
                line = args_getstr(args, "output", 0, "all_sigs.list");
            n = strlen(line);
            i = n;
            tmpstr = malloc(i + 5);
            i--;
            while (i > 0 && line[i]!='.' && line[i]!='/' && line[i]!='\\') {
                i--;
            }
            if (i == 0)i = n;
            memcpy(tmpstr, line, i);
            if (strcmp(args_getstr(args, "keepsigs", 0, "continue"), "text") == 0)
                memcpy(tmpstr + i, ".txt", 5);
            else
                memcpy(tmpstr + i, ".sdb", 5);
            if (PFopen(&f, tmpstr, "wb")) {
                args_report_warning(NULL, "Signature database <%s> could not be created\n", tmpstr);
            }
            else {
                if (strcmp(args_getstr(args, "keepsigs", 0, "continue"), "stop") == 0) {
                    PFclose(f);
                    allsigs = genosig_multiplex(result, nfiles);
                    genosig_export(allsigs, tmpstr);
                    result = genosig_demux(allsigs, &nfiles);
                }
                else if (strcmp(args_getstr(args, "keepsigs", 0, "continue"), "average") == 0) {
                    PFclose(f);
                    allsigs = genosig_multiplex(result, nfiles);
                    sigavg = genosig_avgsig(allsigs, 0);
                    genosig_linkname(sigavg, tmpstr);
                    genosig_export(sigavg, tmpstr);
                    if (sigavg != allsigs) genosig_free(sigavg);
                    result = genosig_demux(allsigs, &nfiles);
                }
                else {
                    for (i = 0;i < nfiles;i++) {
                        if (result[i])
                            genosig_savetxt(result[i], f);
                    }
                    PFclose(f);
                }
            }
        }
        /* cleanup */
        if (taxtree) taxtree_free(taxtree);
        if (refgenes)genosig_free(refgenes);
        refgenes = NULL;
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
    if (strcmp(str, "Same_Subgroup") == 0)return 7;
    if (strcmp(str, "Same_Type") == 0)return 7;
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
        return "Same_Subgroup";
    case 8:
        return "Replicate";
    default:
        return "?";
    }
}

char* GenDisCal_ordered_output_gen_line(genosig_t** signatures,size_t fileA, size_t fileB, int simlevel) {
    char* filename1;
    char* filename2;
    char* output;
    char* simlevelstr;
    size_t k;
    simlevelstr = taxsimtostr(simlevel);
    if (signatures[fileA] && genosig_info_name(signatures[fileA])) filename1 = genosig_info_name(signatures[fileA]);
    else filename1 = "<File not opened>";
    if (signatures[fileB] && genosig_info_name(signatures[fileB])) filename2 = genosig_info_name(signatures[fileB]);
    else filename2 = "<File not opened>";
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
        int issearch, int sort_by_distance, PF_t* outputf, genosig_t** signatures
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
        args_t* args, genosig_t** signatures, size_t nfiles, size_t sig_len,
        int threadid, int numthreads, size_t* p_n, genosig_df mf, any_t metharg) {
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
            simlevel = (int)genodist_bytaxonomy(signatures[ifl_added], signatures[nfiles], metharg);
            if (onlylevel < -1 || simlevel == onlylevel) {
                if (signatures[ifl_added] && signatures[nfiles])
                    dist = mf(signatures[ifl_added], signatures[nfiles], metharg) ;
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
                simlevel = (int)genodist_bytaxonomy(signatures[ifl_added], signatures[ifl_c], metharg);
                if (onlylevel < -1 || simlevel == onlylevel) {
                    if (signatures[ifl_added] && signatures[ifl_c])
                        dist = mf(signatures[ifl_added],signatures[ifl_c], metharg);
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

int GenDisCal_perform_comparisons(args_t* args, genosig_t** signatures, size_t n_search, size_t nfiles, size_t sig_len, int numthreads) {
    genosig_df mf;
    any_t metharg;
    int nocalc;
    int outputmode;
    double binwidth;
    double minval;
    double maxval;
    int onlylevel;
    int issearch;
    int sort_by_distance, do_clustering;
    double* ocmp_results;
    int* ocmp_types;
    datatable_t* otable;
    datatable_t** otable_sub;
    PF_t* outputf;
    char* tmpstr1;
    char* tmpstr2;
    size_t tmpstrlen;
    PF_t* tree_f;
    size_t i, n;
    size_t* p_n;
    size_t nfsquare;
    size_t countsperproc;
    size_t orderedamount;
    size_t nclusters;
    time_t step_start;
    time_t step_duration;
    int hours, minutes, seconds;
    /* parse relevant arguments */
    do_clustering = 0;
    nocalc = set_method_function(args, &mf, &metharg);
    args_report_info(NULL, "Method function selected successfully\n", numthreads);

    outputmode = OUTPUTMODE_NORMAL;
    ocmp_results = NULL;
    ocmp_types = NULL;
    sort_by_distance = args_ispresent(args, "sort_by_distance");
    do_clustering = args_ispresent(args, "clustering");
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
        nfiles-= n_search;
        args_report_info(NULL, "Search mode is enabled for %d signatures\n", (int)n_search);
        orderedamount = nfiles * n_search;
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
            if (signatures[i] && genosig_info_name(signatures[i])) {
                datatable_setcolname(otable, i, genosig_info_name(signatures[i]));
                datatable_setrowname(otable, i, genosig_info_name(signatures[i]));
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
        step_start = time(NULL);
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
               if (outputmode == OUTPUTMODE_MATRIX) {
                    datatable_set(otable, ifl_added, ifl_added, 0.0);
                }
                if (issearch) {
                    for (ifl_c = 0;ifl_c < n_search;ifl_c++) {
                        if (outputmode == OUTPUTMODE_ORDERED) {
                            ocmp_types[ifl_added + ifl_c*nfiles] = -3;
                        }
                        simlevel = (int)genodist_bytaxonomy(signatures[ifl_added], signatures[nfiles + ifl_c], metharg);
                        if (onlylevel < -1 || simlevel == onlylevel) {
                            if (signatures[ifl_added] && signatures[nfiles + ifl_c])
                                dist = mf(signatures[ifl_added], signatures[nfiles+ifl_c], metharg);
                            else
                                dist = 1.0;
                            if (dist >= minval && dist <= maxval) {
                                if (outputmode == OUTPUTMODE_NORMAL) {
                                    if (signatures[ifl_added] && genosig_info_name(signatures[ifl_added])) fn1 = genosig_info_name(signatures[ifl_added]);
                                    else fn1 = "<File could not be opened>";
                                    if (signatures[nfiles] && genosig_info_name(signatures[nfiles + ifl_c])) fn2 = genosig_info_name(signatures[nfiles + ifl_c]);
                                    else fn2 = "<File could not be opened>";
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
                                    datatable_set(otable, ifl_c, ifl_added, dist);
                                }
#pragma omp atomic
                                (*p_n)++;
                                if (nfiles*n_search < 100 || (*p_n) % (nfiles*n_search / 100) == 0) {
#pragma omp critical
                                    args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) files compared\r", (*p_n), nfiles, (int)(((double)n * 100) / (double)(nfiles*n_search)));
                                }
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
                        simlevel = (int)genodist_bytaxonomy(signatures[ifl_added], signatures[ifl_c], metharg);
                        if (onlylevel < -1 || simlevel == onlylevel) {
                            if (signatures[ifl_added] && signatures[ifl_c])
                                dist = mf(signatures[ifl_added], signatures[ifl_c], metharg);
                            else
                                dist = 1.0;
                            #pragma omp atomic
                            (*p_n)++;
                            if (nfsquare < 100 || (*p_n) % (nfsquare / 100) == 0) {
                                #pragma omp critical
                                args_report_progress(NULL, _LLD_ "/" _LLD_ " (%d%%) comparisons performed\r", (*p_n), nfsquare, (int)(((double)n * 100) / (double)nfsquare));
                            }
                            if (dist >= minval && dist <= maxval) {
                                if (outputmode == OUTPUTMODE_NORMAL) {
                                    if (signatures[ifl_added] && genosig_info_name(signatures[ifl_added])) fn1 = genosig_info_name(signatures[ifl_added]);
                                    else fn1 = "<File could not opened>";
                                    if (signatures[ifl_c] && genosig_info_name(signatures[ifl_c])) fn2 = genosig_info_name(signatures[ifl_c]);
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
        step_duration = time(NULL) - step_start;
        hours = (int)(step_duration / (60 * 60));
        minutes = (int)((step_duration % (60 * 60)) / (60));
        seconds = (int)((step_duration % (60)));
        args_report_progress(NULL, "A Total of " _LLD_ " comparisons were performed (%d hr %d min %d sec)\n", n, hours, minutes, seconds);
    }
    tree_f = NULL;
    if (outputmode == OUTPUTMODE_MATRIX && (do_clustering || sort_by_distance)) {
        tmpstr1 = args_getstr(args, "output", 0, "tree");
        tmpstrlen = strlen(tmpstr1);
        tmpstr2 = malloc(tmpstrlen + 5);
        memcpy(tmpstr2, tmpstr1, tmpstrlen);
        i = tmpstrlen;
        while (i > 0 && tmpstr2[i - 1] != '.')i--;
        if (i == 0) i = tmpstrlen;
        else i--;
        memcpy(tmpstr2 + i, ".tre", 5);
        PFopen(&tree_f, tmpstr2, "wb");
        free(tmpstr2);
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
            if (sort_by_distance) {
                step_start = time(NULL);
                if (strcmp(args_getstr(args, "clustering", 0, "single"), "single") == 0) {
                    args_report_progress(NULL, "Sorting the output matrix using single linkage...\n");
                    nclusters = datatable_dmsingleclust(otable, tree_f);
                    do_clustering = 0;
                    PFclose(tree_f);
                }
                else if (strcmp(args_getstr(args, "clustering", 0, "sort"), "ward") == 0) {
                    args_report_warning(NULL, "Ward clustering destroys the distance matrix, and as such the output sorted using single linkage.\n");
                    args_report_progress(NULL, "Sorting the output matrix using single linkage...\n");
                    nclusters = datatable_dmsingleclust(otable, NULL);
                }
                else {
                    args_report_progress(NULL, "Sorting the output matrix using a basic sorting approach...\n");
                    nclusters = datatable_dmrclust(otable, 0.5);
                    do_clustering = 0;
                    PFclose(tree_f);
                }
                step_duration = time(NULL) - step_start;
                hours = (int)(step_duration / (60 * 60));
                minutes = (int)((step_duration % (60 * 60)) / (60));
                seconds = (int)((step_duration % (60)));
                args_report_progress(NULL, "Output matrix seems to have about " _LLD_ " clusters. (%d hr %d min %d sec)\n",
                    nclusters, hours, minutes, seconds);
            }
            datatable_write(otable, outputf, ",", DATATABLE_WRITECOLNAMES | DATATABLE_WRITEROWNAMES);
        }
    }
    if (!nocalc && outputmode == OUTPUTMODE_ORDERED) {
        GenDisCal_print_ordered_output(orderedamount, ocmp_types, ocmp_results, nfiles, issearch, sort_by_distance, outputf, signatures);
    }
    if (!nocalc && outputmode == OUTPUTMODE_MATRIX && do_clustering) {

        args_report_progress(NULL, "Printing dendrogram based on the output matrix...\n");
        step_start = time(NULL);
        if (strcmp(args_getstr(args, "clustering", 0, "ward"), "ward") == 0)
            nclusters = datatable_dmwardclust(otable, 0.01, tree_f);
        if (strcmp(args_getstr(args, "clustering", 0, "ward"), "single") == 0)
            nclusters = datatable_dmsingleclust(otable, tree_f);
        PFclose(tree_f);
        step_duration = time(NULL) - step_start;
        hours = (int)(step_duration / (60 * 60));
        minutes = (int)((step_duration % (60 * 60)) / (60));
        seconds = (int)((step_duration % (60)));
        args_report_progress(NULL, "Output seems to have about " _LLD_ " clusters. (%d hr %d min %d sec)\n", nclusters, hours, minutes, seconds);
    }
    /* cleanup */
    if (otable)datatable_free(otable);
    return 0;
}

int GenDisCal(args_t* args) {
    char** filelist;
    genosig_t** siglist;
    size_t i;
    size_t numfiles, nfilesmod;
    size_t siglen;
    size_t nsearch;
    int numthreads;
    filelist = GenDisCal_get_file_list(args, &numfiles);
    if (args_ispresent(args, "license")) {
        fprintf(stdout, 
                "GenDisCal\n"
                "-------------------------------------------------------------------------------\n"
                "MIT License\n"
                "\n"
                "Copyright (c) 2020 Gleb Goussarov\n"
                "\n"
                "Permission is hereby granted, free of charge, to any person obtaining a copy\n"
                "of this software and associated documentation files (the \"Software\"), to deal\n"
                "in the Software without restriction, including without limitation the rights\n"
                "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
                "copies of the Software, and to permit persons to whom the Software is\n"
                "furnished to do so, subject to the following conditions:\n"
                "\n"
                "The above copyright notice and this permission notice shall be included in all\n"
                "copies or substantial portions of the Software.\n"
                "\n"
                "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
                "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
                "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
                "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
                "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
                "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
                "SOFTWARE.\n"
                "------------------------------------------------------------------------------\n");
    }
    else if (args_ispresent(args, "version")) {
        fprintf(stdout, "%s\n", VERSION_NAME);
    }
    else if (args_ispresent(args, "changes")) {
        fprintf(stdout, "%s\n", CHANGES);
    }
    else if (numfiles < 1) {
        args_report_error(NULL, "Number of files (" _LLD_ ") to be loaded is below 2! Aborting.\n", numfiles);
    }
    else {
        numthreads = args_getint(args, "nprocs", 0, autothreads());
        args_report_info(NULL, "Using up to %d threads\n", numthreads);
        nfilesmod = numfiles;
        siglist = GenDisCal_get_signatures(args, filelist, &nfilesmod, &siglen, numthreads, &nsearch);
        if (nfilesmod < 1) {
            args_report_error(NULL, "Number of files (" _LLD_ ") to be loaded is below 2! Aborting.\n", numfiles);
        }
        else {
            if (!args_ispresent(args, "keepsigs") || strcmp(args_getstr(args, "keepsigs", 0, "continue"), "continue") == 0) {
                GenDisCal_perform_comparisons(args, siglist, nsearch, nfilesmod, siglen, numthreads);
            }

            for (i = 0;i < nfilesmod;i++) {
                if (siglist[i])genosig_free(siglist[i]);
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
