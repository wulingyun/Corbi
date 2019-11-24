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

double cdf_kde(double *x, int n, double v, double bw, double med)
{
  double p, s;
  p = 0;
  if (v <= med)
  {
    for (int i = 0; i < n; i++)
    {
      s = pnorm(v - x[i], 0, bw, 1, 0);
      if (s < 1e-8) break;
      p += s;
    }
    p = p / n;
  }
  else
  {
    for (int i = n-1; i >= 0; i--)
    {
      s = pnorm(v - x[i], 0, bw, 0, 0);
      if (s < 1e-8) break;
      p += s;
    }
    p = (n - p) / n;
  }
  return(p);
}

double quantile_kde(double *x, int n, double p, bool sorted)
{
  if (!sorted) R_rsort(x, n);
  double bw, med, v0, v1, v2, p0, p1, p2;
  bw = bw_nrd0(x, n, true);
  med = quantile(x, n, 0.5, true);
  v0 = quantile(x, n, p, true);
  p0 = cdf_kde(x, n, v0, bw, med);
  if (p0 < p)
  {
    v1 = x[n-1];
    p1 = cdf_kde(x, n, v1, bw, med);
    while (p1 < p)
    {
      v1 = v1 + 1;
      p1 = cdf_kde(x, n, v1, bw, med);
    }
  }
  else
  {
    v1 = v0;
    p1 = p0;
    v0 = x[0];
    p0 = cdf_kde(x, n, v0, bw, med);
    while (p0 > p)
    {
      v0 = v0 - 1;
      p0 = cdf_kde(x, n, v0, bw, med);
    }
  }
  while ((v1 - v0) > 1e-6)
  {
    v2 = (v0 + v1) / 2;
    p2 = cdf_kde(x, n, v2, bw, med);
    if (p2 < p)
    {
      v0 = v2;
      p0 = p2;
    }
    else
    {
      v1 = v2;
      p1 = p2;
    }
  }
  return(v0);
}

