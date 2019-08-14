#include "GenDisCal_distances.h"
#include "vecops.h"
#include "suffixtree.h"
#include <string.h>

/* helper functions */
double C0_step_function(double x, double tmin, double tmax) {
    if (x < tmin)return 0.0;
    else if (x > tmax)return 1.0;
    else if (tmax == tmin)return 0.5;
    else return (x - tmin) / (tmax - tmin);
}
double C0_sigmoid(double x, double t, double factor) {
    return (1.0 / (1.0 + pow(x / t, factor)));
}
double nucfreq(nucseq* src, int nuc) {
    size_t numnuc = 0;
    size_t numfound = 0;
    size_t i_ = 0;
    for (i_ = 0;i_ < src->len;i_++)
    {
        if (src->seq[i_] == nuc)
            numfound++;
        if (nuc == nucN || src->seq[i_] != nucN)
            numnuc++;
    }
    return ((double)numfound) / ((double)numnuc);
}
#ifndef PI
#define PI (atan(1)*4)
#endif
double gaussianbell(double mu, double sigma, double x) {
    double normpos;
    double var;
    normpos = x - mu;
    var = sigma*sigma;
    return exp(-(normpos*normpos)/(2*var))/sqrt(2*PI*var);
}
double cut_gaussianbell(double mu, double sigma, double x, double outlier_range) {
    double gaussvalue;
    double outliervalue;
    gaussvalue = gaussianbell(mu, sigma, x);
    outliervalue = gaussianbell(mu, sigma, mu + outlier_range*sigma);
    if (outliervalue > gaussvalue) return 0;
    else return gaussvalue - outliervalue;
}
// bases
size_t gensz(nucseq* sequence, int n, double** result) {
    *result = malloc(sizeof(double));
    (*result)[0] = (double) (sequence->len);
    return sizeof(double) * 1;
}
size_t freqn(nucseq* sequence, int n, double** result) { // n-mer frequency signature
    int32_t* counts;
    counts = oligocount(sequence, n);
    *result = freqsig(counts, n);
    free(counts);
    return ((size_t)1 << (2 * n)) * sizeof(double);
}
size_t freqs(nucseq* sequence, int n, double** result) { // n-mer frequency signature
    int32_t* counts;
    counts = oligocount_2strand(sequence, n);
    *result = freqsig(counts, n);
    free(counts);
    return ((size_t)1 << (2 * n)) * sizeof(double);
}
size_t markz1(nucseq* sequence, int n, double** result) {
    int32_t* counts;
    int c;
    size_t comb, nucmask;
    double fs[4];
    counts = oligocount(sequence, n);
    fs[nucA] = nucfreq(sequence, nucA);
    fs[nucC] = nucfreq(sequence, nucC);
    fs[nucG] = nucfreq(sequence, nucG);
    fs[nucT] = nucfreq(sequence, nucT);
    *result = freqsig(counts, n);
    for (comb = 0;comb < (size_t)1 << (2 * n); comb++) {
        for (c = 0; c < n;c++) {
            nucmask = 0b11 << (2 * c);
            (*result)[comb] /= fs[(comb&nucmask) >> (2 * c)];
        }
    }
    free(counts);
    return ((size_t)1 << (2 * n)) * sizeof(double);
}
size_t markz2(nucseq* sequence, int n, double** result) {
    int32_t* counts;
    int c;
    size_t comb, nucmask;
    double gc, fs[4];
    counts = oligocount_2strand(sequence, n);
    gc = nucseq_GC(sequence);
    fs[nucA] = (1.0 - gc) / 2.0;
    fs[nucC] = gc / 2.0;
    fs[nucG] = fs[nucC];
    fs[nucT] = fs[nucA];
    *result = freqsig(counts, n);
    for (comb = 0;comb < (size_t)1 << (2 * n); comb++) {
        for (c = 0; c < n;c++) {
            nucmask = 0b11 << (2 * c);
            (*result)[comb] /= fs[(comb&nucmask) >> (2 * c)];
        }
    }
    free(counts);
    return ((size_t)1 << (2 * n)) * sizeof(double);
}
size_t TETRA(nucseq* sequence, int unused, double** result) {
    int32_t* counts1;
    int32_t* counts2;
    size_t i, numel;
    nucseq revcomp = EMPTYSEQ;
    counts1 = oligocount(sequence, 4);
    nucseqrevcomp(sequence, &revcomp);
    counts2 = oligocount(&revcomp, 4);
    numel = 256;
    for (i = 0;i < numel;i++)
        counts1[i] += counts2[i];
    *result = TETRAsig(counts1);
    free(counts1);
    free(counts2);
    clear_nucseq(&revcomp);
    return 256 * sizeof(double);
}
size_t karln(nucseq* sequence, int n, double** result) { // n-mer karlin signature
    int32_t* counts;
    counts = oligocount(sequence, n);
    *result = karlinsig(counts, n);
    free(counts);
    return ((size_t)1 << (2 * n)) * sizeof(double);
}
size_t karsn(nucseq* sequence, int n, double** result) { // n-mer karlin* signature
    int32_t* counts;
    size_t numel;
    nucseq revcomp = EMPTYSEQ;
    counts = oligocount_2strand(sequence, n);
    *result = karlinsig(counts, n);
    free(counts);
    numel = (size_t)1 << (n * 2);
    return numel * sizeof(double);
}
size_t markn1(nucseq* sequence, int n, double** result) {
    int32_t* counts;
    size_t numel;
    nucseq revcomp = EMPTYSEQ;
    counts = oligocount(sequence, n);
    *result = fullmarkovsig(counts, n);
    free(counts);
    numel = (size_t)1 << (n * 2);
    return numel * sizeof(double);
}
size_t markn2(nucseq* sequence, int n, double** result){
    int32_t* counts;
    size_t numel;
    nucseq revcomp = EMPTYSEQ;
    counts = oligocount_2strand(sequence, n);
    *result = fullmarkovsig(counts, n);
    free(counts);
    numel = (size_t)1 << (n * 2);
    return numel * sizeof(double);
}
size_t multikarl(nucseq* sequence, int n, double** result) {
    *result = multikarlsig(&sequence, 1, n, 20000, 10000);
    return (1LL << (2 * n))*(1LL << (2 * n)) * sizeof(double);
}
size_t multifreq(nucseq* sequence, int n, double** result) {
    *result = multifreqsig(&sequence, 1, n, 3000, "ATG");
    return (1LL << (2 * n))*(1LL << (2 * n)) * sizeof(double);
}
size_t minhashsig(nucseq* sequence, int n, double** result) {
    *result = (double*) minhash_Msig(&sequence, 1, n, SIGLEN_MINHASH, 200000);
    return SIGLEN_MINHASH*sizeof(int64_t);
}

size_t combinedsig(nucseq* sequence, int n, double ** result)
{
    double *K4;
    double *K6;
    double *minhash;
    double *lensig;
    size_t sigsize4, sigsize6, sigsizem, sigsizel;
    char* tmpres;
    sigsize4 = karsn(sequence, 4, &K4);
    sigsize6 = karsn(sequence, 6, &K6);
    sigsizem = minhashsig(sequence, 31, &minhash);
    sigsizel = gensz(sequence, 0, &lensig);
    tmpres = malloc(sigsize4 + sigsize6 + sigsizem);
    memcpy(tmpres, K4, sigsize4);
    memcpy(tmpres + sigsize4, K6, sigsize6);
    memcpy(tmpres + sigsize4 + sigsize6, minhash, sigsizem);
    memcpy(tmpres + sigsize4 + sigsize6 + sigsizem, minhash, sigsizem);
    *result = (double*) tmpres;
    free(K4);
    free(K6);
    free(minhash);
    free(lensig);
    return sigsize4 + sigsize6 + sigsizem + sigsizel;
}

size_t SSUseq(nucseq* sequence, int unused, double ** result)
{
    static char* consensus = 
        "TAATTGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGTATTGAAAGG"
        "TGCTTGCACCTGGACGAGTGGCGGACGGGTGAGTAACACGTGGGAACCTGCCCTTAAGTGGGGGATAACATTTGGAAA"
        "CAGATGCTAATACCGCATAAAACCGCATGGTTAAAGGCTGAAAGTGGCGTGAGCTATCGCTTTTGGATGGGCCCGCGT"
        "CGGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCGTAGCCGGTCTGAGAGGATGATCGGCCACACT"
        "GGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGC"
        "AATGCCGCGTGAGTGAAGAAGGCCTTCGGGTTGTAAAGCTCTTTTGTTGGGGAAGAAAGGTCGGCAGGTAACTGTTGT"
        "CGGCGTGACGGTACCCAACGAAAAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTT"
        "ATCCGGAATTATTGGGCGTAAAGCGAGCGCAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCTCGGCTTAACCGGGGAA"
        "GTGCATTGGAAACTGGGAAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATA"
        "TGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"
        "ATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGCTAGGTGTTGGGGGGTTTCCGCCCTTCAGTGCCGCAGCT"
        "AACGCATTAAGCATTCCGCCTGGGGAGTACGGCCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCG"
        "GTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCTTTTGAACACCTTAGAGATAA"
        "GGTTTTCCCTTCGGGGACAAAATGACAGGTGCTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCC"
        "CGCAACGAGCGCAACCCTTATCATTAGTTGCCAGCATTAAGTTGGGCACTCTAGTGAGACTGCCGGTGACAAACCGGA"
        "GGAAGGTGGGGATGACGTCAAATCATCATGGCCCTTATGACCTGGGCTACACACGTGCTACAATGGATGGTACAAAGA"
        "GTTGCGAGACCGCGAGGTCAAGCTAATCTCTTAAAGCCATTCTCAGTTCGGATTGTAGTCTGCAACTCGACTACATGA"
        "AGTCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACA"
        "CCATGGGAGTTTGTAACACCCGAAGTCGGTGGCCTAACCTTAGGGAGGGAGCCGACTAAGGTGGGACAGATGACTGGG"
        "GTGAAGTCGTAACAAGGTAGCCGTAGGGGAACCTGCGGCTGGATCACCTCCTTTCT";
    static size_t consensuslen = 0;
    static suftree_t* consensustree = NULL;
    char* marks;
    char** detectedSSU;
    char** detectedrcSSU;
    double* alignment_score;
    size_t* SSUlen;
    size_t* rcSSUlen;
    size_t i, nSSU, nrcSSU, maxlenid, resultlen;
    nucseq tmpseq = EMPTYSEQ;
    
    resultlen = 0;
    if (!consensuslen)
        consensuslen = strlen(consensus);
    if (!consensustree) {
        nucseq_from_string(&tmpseq, consensus);
        consensustree = suftree_from(tmpseq.seq, SUFTREE_DEFAULT, consensuslen);
        clear_nucseq(&tmpseq);
    }
    if (!sequence) {
        suftree_free(consensustree);
        consensustree = NULL;
    }
    else {
        marks = suftree_approximatesearch(sequence->seq, sequence->len, consensustree, 14, 1100, 400);
        detectedSSU = marks_split(sequence->seq, marks, sequence->len, &SSUlen, &nSSU, KEEP_MARKED);
        nucseqrevcomp(sequence, &tmpseq);
        free(marks);
        marks = suftree_approximatesearch(tmpseq.seq, tmpseq.len, consensustree, 14, 1100, 400);
        detectedrcSSU = marks_split(tmpseq.seq, marks, tmpseq.len, &rcSSUlen, &nrcSSU, KEEP_MARKED);
        free(marks);
        if (nSSU == 0 && nrcSSU > 0) {
            nSSU = nrcSSU;
            detectedSSU = detectedrcSSU;
            SSUlen = rcSSUlen;
        }
        else if (nSSU > 0 && nrcSSU > 0) {
            detectedSSU = (char**)realloc(detectedSSU, sizeof(char*)*(nSSU + nrcSSU));
            SSUlen = (size_t*)realloc(SSUlen, sizeof(size_t)*(nSSU + nrcSSU));
            memcpy(detectedSSU + nSSU, detectedrcSSU, sizeof(char*)*nrcSSU);
            memcpy(SSUlen + nSSU, rcSSUlen , sizeof(size_t)*nrcSSU);
            free(detectedrcSSU);
            free(rcSSUlen);
            nSSU += nrcSSU;
        }
        if (nSSU > 0) {
            alignment_score = (double*)malloc(sizeof(double)*nSSU);
            alignment_score[0] = 0.0;
            maxlenid = 0;
            for (i = 0;i < nSSU;i++) {
                alignment_score[i] = suftree_roughalign(detectedSSU[i], SSUlen[i], consensustree, 8, 1, 1);
                if (alignment_score[i] < alignment_score[maxlenid])maxlenid = i;
            }
            *result = (double*)suftree_from(detectedSSU[maxlenid], SUFTREE_DEFAULT, SSUlen[maxlenid]);
            resultlen = SSUlen[maxlenid];
            for (i = 0;i < nSSU;i++) {
                free(detectedSSU[i]);
            }
            free(alignment_score);
            free(detectedSSU);
            free(SSUlen);
        }
    }
    return resultlen;
}

// methods
double reldist(double* sig1, double* sig2, size_t veclen, double unused) {
    double* tmp;
    double result;
    size_t i;
    tmp = malloc(sizeof(double)*veclen);
    for (i = 0;i < veclen;i++) {
        tmp[i] = (sig1[i] - sig2[i])*2.0 / (sig1[i] + sig2[i]);
    }
    vec_abs(tmp,veclen);
    result = vec_avg(tmp, veclen);
    free(tmp);
    return result;
}
double ED(double* sig1, double* sig2, size_t veclen, double unused) {
    double* tmp;
    double result;
    tmp = malloc(sizeof(double)*veclen);
    vec_zero(tmp, veclen);
    vec_add(tmp, sig1, veclen);
    vec_subtract(tmp, sig2, veclen);
    result = vec_norm(tmp, veclen);
    free(tmp);
    return result;
}
double AMD(double* sig1, double* sig2, size_t veclen, double unused) {
    double* tmp;
    double result;
    tmp = malloc(sizeof(double)*veclen);
    vec_zero(tmp, veclen);
    vec_add(tmp, sig1, veclen);
    vec_subtract(tmp, sig2, veclen);
    vec_abs(tmp, veclen);
    result = vec_avg(tmp, veclen);
    free(tmp);
    return result;
}
double SVC(double* sig1, double* sig2, size_t veclen, double threshold) { // Count similar values
    double* tmp;
    double result;
    double t1;
    double t2;
    size_t i;
    // before we begin, we need to find reasonable values for t1 and t2
    // for the sake of simplicity, we take t1 = 0.0 and t2 = threshold
    t1 = 0.0;
    t2 = threshold;
    result = 0;
    tmp = malloc(sizeof(double)*veclen);
    vec_zero(tmp, veclen);
    vec_add(tmp, sig1, veclen);
    vec_subtract(tmp, sig2, veclen);
    vec_abs(tmp, veclen);
    for (i = 0;i < veclen;i++)
        result += C0_step_function(tmp[i], t1, t2);
    result /= (double)veclen;
    free(tmp);
    return result;
}
double SSVC(double* sig1, double* sig2, size_t veclen, double threshold) { // Count similar values
    double* tmp;
    double result;
    double t1;
    double t2;
    size_t i;
    // before we begin, we need to find reasonable values for t1 and t2
    // for the sake of simplicity, we take t1 = 0.0 and t2 = threshold
    t1 = 0.0;
    t2 = threshold;
    result = 0;
    tmp = malloc(sizeof(double)*veclen);
    vec_zero(tmp, veclen);
    vec_add(tmp, sig1, veclen);
    vec_subtract(tmp, sig2, veclen);
    vec_dot(tmp, tmp, veclen);
    for (i = 0;i < veclen;i++)
        result += C0_step_function(tmp[i], t1, t2);
    result /= (double)veclen;
    free(tmp);
    return result;
}
double same_species(double* sig1, double* sig2, size_t veclen, double code) {
    double value;
    if (code == 1.0) {
        value = SVC(sig1, sig2, veclen, 0.02);
        if (value < 0.10)return 0.99;
        if (value > 0.30)return 0.01;
        return 0.5;
    }
    return 0.0;
}
double ESVC(double* sig1, double* sig2, size_t veclen, double threshold) {
    double* tmp;
    double result;
    double t1;
    double t2;
    double compdiff;
    size_t i;
    // before we begin, we need to find reasonable values for t1 and t2
    // for the sake of simplicity, we take t1 = 0.0 and t2 = threshold
    t1 = 0.0;
    t2 = threshold;
    result = 0;
    tmp = malloc(sizeof(double)*veclen);
    vec_zero(tmp, veclen);
    vec_add(tmp, sig1, veclen);
    vec_subtract(tmp, sig2, veclen);
    for (i = 0;i < veclen;i++) {
        compdiff = C0_step_function(tmp[i], t1, t2);
        result += compdiff*compdiff;
    }
    result /= (double)(veclen)*(double)(veclen);
    result = sqrt(result);
    free(tmp);
    return result;
}
double DVC(double* sig1, double* sig2, size_t veclen, double threshold) { // Count similar values
    double* tmp;
    double result;
    size_t i;
    // before we begin, we need to find reasonable values for t1 and t2
    // for the sake of simplicity, we take t1 = 0.0 and t2 = threshold
    result = 0;
    tmp = malloc(sizeof(double)*veclen);
    vec_zero(tmp, veclen);
    vec_add(tmp, sig1, veclen);
    vec_subtract(tmp, sig2, veclen);
    for (i = 0;i < veclen;i++)
        result += C0_sigmoid(tmp[i], threshold, 4);
    result /= (double)veclen;
    free(tmp);
    return result;
}
double MD(double* sig1, double* sig2, size_t veclen, double unused) {
    double* tmp;
    double result;
    tmp = malloc(sizeof(double)*veclen);
    vec_zero(tmp, veclen);
    vec_add(tmp, sig1, veclen);
    vec_subtract(tmp, sig2, veclen);
    vec_abs(tmp, veclen);
    result = vec_max(tmp, veclen);
    free(tmp);
    return result;
} // maximal difference
/* pearson correclation coefficient */
double corr(double* sig1, double* sig2, size_t veclen, double unused) {
    double* tmp1;
    double* tmp2;
    double* tmp3;
    double result;
    double avg1;
    double avg2;
    double num;
    double var1;
    double var2;
    double denom;
    tmp1 = malloc(sizeof(double)*veclen);
    tmp2 = malloc(sizeof(double)*veclen);
    tmp3 = malloc(sizeof(double)*veclen);

    vec_zero(tmp1, veclen);
    vec_add(tmp1, sig1, veclen);
    avg1 = vec_avg(tmp1, veclen);
    vec_subtract_all(tmp1, avg1, veclen);

    vec_zero(tmp2, veclen);
    vec_add(tmp2, sig2, veclen);
    avg2 = vec_avg(tmp2, veclen);
    vec_subtract_all(tmp2, avg2, veclen);

    vec_zero(tmp3, veclen);
    vec_add(tmp3, tmp1, veclen);
    vec_dot(tmp3, tmp2, veclen);
    num = vec_sum(tmp3, veclen);

    vec_dot(tmp1, tmp1, veclen);
    vec_dot(tmp2, tmp2, veclen);
    if (veclen > 1) {
        var1 = vec_sum(tmp1, veclen);
        var2 = vec_sum(tmp2, veclen);
    }
    else {
        var1 = var2 = 1;
    }
    denom = sqrt(var1*var2);

    result = 0.5 - (0.5*num / denom);
    free(tmp1);
    free(tmp2);
    free(tmp3);
    return result;
}
double sqrtcorr(double* sig1, double* sig2, size_t veclen, double unused) {
    return sqrt(corr(sig1, sig2, veclen, unused));
}
double SVcorr(double* sig1, double* sig2, size_t veclen, double threshold) { // Correlate similar values
    double* tmp1;
    double* tmp2;
    double* tmp3;
    double result;
    double avg1;
    double avg2;
    double num;
    double var1;
    double var2;
    double denom;
    double t1;
    double t2;
    double dt;
    double diff, midpointx2;
    size_t i;
    t1 = threshold / 2;
    t2 = threshold;
    dt = t2 - t1;
    tmp1 = malloc(sizeof(double)*veclen);
    tmp2 = malloc(sizeof(double)*veclen);
    tmp3 = malloc(sizeof(double)*veclen);

    vec_zero(tmp1, veclen);
    vec_add(tmp1, sig1, veclen);
    avg1 = vec_avg(tmp1, veclen);
    vec_subtract_all(tmp1, avg1, veclen);

    vec_zero(tmp2, veclen);
    vec_add(tmp2, sig2, veclen);
    avg2 = vec_avg(tmp2, veclen);
    vec_subtract_all(tmp2, avg2, veclen);

    for (i = 0;i < veclen;i++) {
        diff = tmp1[i] - tmp2[i];
        midpointx2 = tmp1[i] + tmp2[i];
        if (diff > 0.0) {
            diff = C0_step_function(diff, t1, t2)*dt;
            tmp1[i] = (midpointx2 + diff) / 2;
            tmp2[i] = (midpointx2 - diff) / 2;
        }
        else if (diff < 0.0) {
            diff = -C0_step_function(diff, t1, t2)*dt;
            tmp1[i] = (midpointx2 - diff) / 2;
            tmp2[i] = (midpointx2 + diff) / 2;
        }
    }

    vec_zero(tmp3, veclen);
    vec_add(tmp3, tmp1, veclen);
    vec_dot(tmp3, tmp2, veclen);
    num = vec_sum(tmp3, veclen);


    vec_dot(tmp1, tmp1, veclen);
    vec_dot(tmp2, tmp2, veclen);
    if (veclen > 1) {
        var1 = vec_sum(tmp1, veclen);
        var2 = vec_sum(tmp2, veclen);
    }
    else {
        var1 = var2 = 1;
    }
    denom = sqrt(var1*var2);

    result = 0.5 - (0.5*num / denom);
    free(tmp1);
    free(tmp2);
    free(tmp3);
    return result;
}
double multiminSVC(double* sig1, double* sig2, size_t veclen, double threshold) {
    double result;
    double tmpresult;
    size_t subveclen;
    size_t i;
    subveclen = (size_t)round(sqrt((double)veclen));
    result = SVC(sig1, sig2, subveclen, threshold);
    for (i = 1;i < subveclen;i++) {
        tmpresult = SVC(sig1 + i*subveclen, sig2 + i*subveclen, subveclen, threshold);
        if (tmpresult < result)result = tmpresult;
    }
    return result;
}
double SNPest(double* freq1, double* freq2, size_t kmerlen, size_t seqlen) {
    size_t i;
    double result = 0;
    double delta;
    size_t veclen = 1LL << (2 * kmerlen);
    double errorrate = (double)(2 * kmerlen);
    for (i = 0;i < veclen;i++) {
        delta = freq1[i] - freq2[i];
        if (delta < 0)
            delta = -delta;
        result += delta;
    }
    result = result/errorrate;
    return result;
}
double multifreqdist(double* sig1, double* sig2, size_t veclen, double threshold) {
    double result;
    double tmpresult;
    double* dists;
    size_t subveclen;
    size_t i,j;
    size_t kmerlen;
    subveclen = (size_t)round(sqrt((double)veclen));
    dists = calloc(subveclen, sizeof(double));
    kmerlen = (size_t)(log2((double)subveclen)/2);
    j = 0;
    for (i = 1;i < subveclen;i++) {
        tmpresult = SNPest(sig1 + i*subveclen, sig2 + i*subveclen, kmerlen, 3000);
        if (tmpresult <= threshold) {
            dists[j] = tmpresult;
            j++;
        }
    }
    result = vec_avg(dists,j);
    free(dists);
    return result;
}
static inline double resemblance2ANI(double v, double k, double genomesize) {
    if (v == 0.0)return 0;
    return 1.0+log(v*2/(1+v))/k;
}
double approxANI(double* sig1, double* sig2, size_t veclen, double k) {
    int64_t* s1;
    int64_t* s2;
    size_t i,j;
    size_t AinterB;
    size_t AunionB;
    double resemblance;
    s1 = (int64_t*)sig1;
    s2 = (int64_t*)sig2;
    i = j = 0;
    AinterB = AunionB = 0;
    while (i < veclen && j < veclen) {
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
    return 1-resemblance2ANI(resemblance, k, (double)((s1[0]+s2[0])/2));
}

double combinedSpecies(double* sig1, double* sig2, size_t veclen, double k)
{
    /* For PaSiT4, the distributions are modelled as normals, with parameters mu and sigma:
         	Species	Genus	Family	Order	Class	Phylum	Kingdom
       MU	0.17	0.75	0.82	0.85	0.87	0.89	0.89
    SIGMA	0.079	0.083	0.068	0.051	0.039	0.036	0.031
       We consider that Same Species, Same Genus or Different Genus are equally likely to occur,
       and approximate "Different Genus" as the "Same Family" distribution in this case
       
       For ANI, the following values were used
         	Species	Different species
       MU	0.02	0.15	0.20
    SIGMA	0.01	0.03	0.10
       The "Different Species" distirbution is a bit difficult to evaluate in this case so
       the value is somewhat arbitrary
    */
    double* K4A;
    double* K6A;
    double* minhashA;
    double* K4B;
    double* K6B;
    double* minhashB;
    double PaSiT4value, PaSiT4conclusion;
    double minhashvalue, minhashconclusion;
    double spweight, gweight, ngweight;
    
    K4A = sig1;
    K4B = sig2;
    K6A = sig1+256;
    K6B = sig2+256;
    minhashA = sig1+4096+256;
    minhashB = sig2+4096+256;

    PaSiT4value = SVC(K4A, K4B, 256, 0.02);
    spweight = cut_gaussianbell(0.17, 0.079, PaSiT4value, 5);
    gweight = cut_gaussianbell(0.74, 0.083, PaSiT4value, 5);
    ngweight = cut_gaussianbell(0.82, 0.068, PaSiT4value, 5);
    PaSiT4conclusion = spweight / (spweight + gweight + ngweight);

    minhashvalue = approxANI(minhashA, minhashB, SIGLEN_MINHASH, k);
    spweight = gaussianbell(0.02, 0.01, minhashvalue);
    gweight = gaussianbell(0.15, 0.025, minhashvalue);
    minhashconclusion = spweight / (spweight + gweight*2);
    /*
    Weights were derived from the Matthew's corelation coeffcient
    */
    return 1 - (0.76*minhashconclusion + 0.74*PaSiT4conclusion)/(0.76+0.74);
}

double localalign(double * sig1, double * sig2, size_t veclen, double threshold)
{
    suftree_t* s1;
    suftree_t* s2;
    char* seq2;
    size_t seq2len;
    s1 = (suftree_t*)sig1;
    s2 = (suftree_t*)sig2;
    seq2 = suftree_getstrptr(s2, &seq2len);
    return suftree_roughalign(seq2, seq2len, s1, 9, 1, 1);
}
