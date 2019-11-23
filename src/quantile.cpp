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

double var(double *x, int n)
{
  double v = 0;
  if (n > 1)
  {
    double m = 0;
    for (int i = 0; i < n; i++)
      m += x[i];
    m /= n;
    for (int i = 0; i < n; i++)
    {
      double s = x[i] - m;
      v += s * s;
    }
    v /= (n-1);
  }
  else
    v = 0;
  return(v);
}

double bw_nrd0(double *x, int n, bool sorted)
{
  if (!sorted) R_rsort(x, n);
  double sd, IQR, s;
  sd = sqrt(var(x, n));
  IQR = quantile(x, n, 0.75, true) - quantile(x, n, 0.25, true);
  s = min(sd, IQR / 1.34);
  if (s == 0) s = sd;
  if (s == 0) s = fabs(x[0]);
  if (s == 0) s = 1;
  s = 0.9 * s * pow(n, -0.2);
  return(s);
}
