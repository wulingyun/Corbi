#include <R.h>
#include <Rdefines.h>

/* median */

int median(int *v, int n);
double median(double *v, int n);

/* quantile */

double quantile(int *x, int n, double p, bool sorted = false);
double quantile(double *x, int n, double p, bool sorted = false);
