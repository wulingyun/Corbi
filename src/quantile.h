#include <R.h>
#include <Rdefines.h>

/* median */

int median(int *v, int n);
double median(double *v, int n);

/* quantile */

double quantile(int *x, int n, double p, bool sorted = false);
double quantile(double *x, int n, double p, bool sorted = false);

double var(double *x, int n);
double bw_nrd0(double *x, int n, bool sorted = false);
double cdf_kde(double *x, int n, double v, double bw, double m);
double quantile_kde(double *x, int n, double p, bool sorted = false);
