#include "Corbi.h"

/* median */

int median(int *x, int n)
{
  R_isort(x, n);
  if (n % 2 == 1)
    return x[n/2+1];
  else
    return (x[n/2] + x[n/2+1])/2;
}

double median(double *x, int n)
{
  R_rsort(x, n);
  if (n % 2 == 1)
    return x[n/2];
  else
    return (x[n/2-1] + x[n/2])/2;
}

double quantile(int *x, int n, double p, bool sorted)
{
  if (!sorted) R_isort(x, n);
  int j;
  double q, m, gamma;
  m = 1-p;
  j = floor(n*p+m);
  gamma = n*p+m-j;
  q = (1-gamma) * x[j-1];
  if (j < n) q += gamma * x[j];
  return (q);
}

double quantile(double *x, int n, double p, bool sorted)
{
  if (!sorted) R_rsort(x, n);
  int j;
  double q, m, gamma;
  m = 1-p;
  j = floor(n*p+m);
  gamma = n*p+m-j;
  q = (1-gamma) * x[j-1];
  if (j < n) q += gamma * x[j];
  return (q);
}
