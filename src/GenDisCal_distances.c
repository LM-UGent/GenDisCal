#include "GenDisCal_distances.h"
#include "vecops.h"

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

// bases
size_t freqn(nucseq* sequence, int n, double** result) { // n-mer frequency signature
    int32_t* counts;
    counts = oligocount(sequence, n);
    *result = freqsig(counts, n);
    free(counts);
    return (size_t)1 << 2 * n;
}
size_t freqs(nucseq* sequence, int n, double** result) { // n-mer frequency signature
    int32_t* counts;
    counts = oligocount_2strand(sequence, n);
    *result = freqsig(counts, n);
    free(counts);
    return (size_t)1 << 2 * n;
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
    return (size_t)1 << (2 * n);
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
    return (size_t)1 << (2 * n);
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
    return 256;
}
size_t karln(nucseq* sequence, int n, double** result) { // n-mer karlin signature
    int32_t* counts;
    counts = oligocount(sequence, n);
    *result = karlinsig(counts, n);
    free(counts);
    return (size_t)1 << 2 * n;
}
size_t karsn(nucseq* sequence, int n, double** result) { // n-mer karlin* signature
    int32_t* counts;
    size_t numel;
    nucseq revcomp = EMPTYSEQ;
    counts = oligocount_2strand(sequence, n);
    *result = karlinsig(counts, n);
    free(counts);
    numel = (size_t)1 << (n * 2);
    return numel;
}
size_t markn1(nucseq* sequence, int n, double** result) {
    /* TODO: not implemented*/
    return 0;
}
size_t markn2(nucseq* sequence, int n, double** result){
    /* TODO: not implemented*/
    return 0;
}
size_t multikarl(nucseq* sequence, int n, double** result) {
    *result = multikarlsig(&sequence, 1, n, 20000, 10000);
    return (1LL << (2 * n))*(1LL << (2 * n));
}
// methods
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