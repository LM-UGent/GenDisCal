#include <math.h>
#include "stats.h"
#include "vecops.h"

/* functions */
double mfunc_beta(double x, double y) {
    return exp(lgamma(x)+lgamma(y)-lgamma(x + y));
}
double mfunc_reg_incbeta(double x, double a, double b) {
    return mfunc_incbeta(x, a, b) / mfunc_beta(a, b);
}
/* "test" functions return the values of the CDF */
double test_snedecor(double x, double df_num, double df_denom) {
    return mfunc_reg_incbeta(df_num*x / (df_num*x + df_denom), df_num / 2, df_denom / 2);
}