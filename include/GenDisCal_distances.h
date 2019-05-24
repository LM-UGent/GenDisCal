#ifndef DEF_GENDISCAL_BM
#define DEF_GENDISCAL_BM

#include "nucseq.h"

typedef size_t(*basis_function)(nucseq*, int, double**);

size_t freqn(nucseq* sequence, int n, double** result);
size_t freqs(nucseq* sequence, int n, double** result);
size_t TETRA(nucseq* sequence, int unused, double** result);
size_t karln(nucseq* sequence, int n, double** result);
size_t karsn(nucseq* sequence, int n, double** result);
size_t markz1(nucseq* sequence, int n, double** result);
size_t markz2(nucseq* sequence, int n, double** result);
size_t markn1(nucseq* sequence, int n, double** result);
size_t markn2(nucseq* sequence, int n, double** result);
size_t multikarl(nucseq* sequence, int n, double** result);
size_t multifreq(nucseq* sequence, int n, double** result);
#define SIGLEN_MINHASH  2000
size_t minhashsig(nucseq* sequence, int n, double** result);

typedef double(*method_function)(double*, double*, size_t, double);

double ED(double* sig1, double* sig2, size_t veclen, double unused);
double AMD(double* sig1, double* sig2, size_t veclen, double unused);
double SVC(double* sig1, double* sig2, size_t veclen, double threshold);
double SSVC(double* sig1, double* sig2, size_t veclen, double threshold);
double same_species(double* sig1, double* sig2, size_t veclen, double code);
double ESVC(double* sig1, double* sig2, size_t veclen, double threshold);
double DVC(double* sig1, double* sig2, size_t veclen, double threshold);
double MD(double* sig1, double* sig2, size_t veclen, double unused);
  /* pearson correclation coefficient */
double corr(double* sig1, double* sig2, size_t veclen, double unused);
double sqrtcorr(double* sig1, double* sig2, size_t veclen, double unused);
double SVcorr(double* sig1, double* sig2, size_t veclen, double threshold);
double multiminSVC(double* sig1, double* sig2, size_t veclen, double threshold);
double approxANI(double* sig1, double* sig2, size_t veclen, double threshold);

#endif

