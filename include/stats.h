#ifndef DEF_STATS
#define DEF_STATS
#include <math.h>
#include "incbeta.h"

/* functions */
#define  mfunc_gamma tgamma
double mfunc_beta(double x, double y);
#define mfunc_incbeta incbeta
double mfunc_reg_incbeta(double x, double a, double b);

double test_snedecor(double x, double df_num, double df_denom);

#endif

