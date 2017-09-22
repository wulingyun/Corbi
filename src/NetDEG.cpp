#include "Corbi.h"
#include <Rmath.h>

SEXP ND_PvalueNetDEG(SEXP _NetDegree, SEXP _nGenes, SEXP _pEdge)
{
  PROTECT(_NetDegree = AS_INTEGER(_NetDegree));
  int *NetDegree = INTEGER_POINTER(_NetDegree);
  int nNetDegree = length(_NetDegree);

  int nGenes = INTEGER_POINTER(AS_INTEGER(_nGenes))[0];
  double pEdge = NUMERIC_POINTER(AS_NUMERIC(_pEdge))[0];

  double *fast_dbinom_k = (double *) R_alloc(nGenes, sizeof(double));
  for (int i = 0; i < nGenes; i++)
  {
    fast_dbinom_k[i] = dbinom(i, nGenes-1, pEdge, 0);
  }
  
  double **fast_dbinom_d = (double **) R_alloc(nGenes, sizeof(double *));
  for (int k = 0; k < nGenes; k++)
  {
    int kd = floor(0.5*k) + 1;
    fast_dbinom_d[k] = (double *) R_alloc(kd, sizeof(double));
    fast_dbinom_d[k][0] = dbinom(0, k, 0.5, 0);
    for (int i = 1; i < kd; i++)
    {
      fast_dbinom_d[k][i] = fast_dbinom_d[k][i-1] + dbinom(i, k, 0.5, 0);
    }
  }

  SEXP _P;
  PROTECT(_P = NEW_NUMERIC(nNetDegree));
  double *P = NUMERIC_POINTER(_P);
  
  for (int i = 0; i < nNetDegree; i++)
  {
    P[i] = 0;
    int d = NetDegree[i];
    for (int k = d; k < nGenes; k++)
    {
      int j = (int) floor(0.5*(k-d));
      P[i] += NetDegree[i] > 0 ? fast_dbinom_k[k] * fast_dbinom_d[k][j] : fast_dbinom_k[k] * 0.5;
    }
  }

  UNPROTECT(2);
  return(_P);
}

SEXP ND_RatioDistribution(SEXP _ExprVal)
{
  PROTECT(_ExprVal = AS_NUMERIC(_ExprVal));
  double *ExprVal = NUMERIC_POINTER(_ExprVal);
  int *dim = INTEGER_POINTER(AS_INTEGER(GET_DIM(_ExprVal)));
  int nGenes = dim[0];
  int nSamples = dim[1];

  SEXP _Median;
  PROTECT(_Median = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_Median, nGenes, nGenes);
  double *Median = NUMERIC_POINTER(_Median);

  SEXP _MAD;
  PROTECT(_MAD = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_MAD, nGenes, nGenes);
  double *MAD = NUMERIC_POINTER(_MAD);

  for (int i = 0; i < nGenes * nGenes; i++)
  {
    Median[i] = 0;
    MAD[i] = 0;
  }

  double *r = (double *) R_alloc(nSamples, sizeof(double));
  double *e, m;
  int n;
  for (int i = 0; i < nGenes-1; i++)
  {
    for (int j = i+1; j < nGenes; j++)
    {
      e = ExprVal;
      n = 0;
      for (int k = 0; k < nSamples; k++)
      {
        if (!ISNA(e[i]) && !ISNA(e[j]) && !ISNAN(e[i]) && !ISNAN(e[j]))
          r[n++] = (e[i] - e[j]) / (e[i] + e[j]);
        e += nGenes;
      }

      m = median(r, n);
      Median[i+nGenes*j] = m;
      Median[j+nGenes*i] = -m;

      for (int k = 0; k < n; k++)
        r[k] = fabs(r[k] - m);

      m = median(r, n) * 1.4826;
      MAD[i+nGenes*j] = m;
      MAD[j+nGenes*i] = m;
    }
  }

  SEXP _ratio;
  PROTECT(_ratio = NEW_LIST(2));
  SetListElement(_ratio, 0, "median", _Median);
  SetListElement(_ratio, 1, "mad", _MAD);

  UNPROTECT(4);
  return(_ratio);
}

SEXP ND_RatioNet(SEXP _RatioMedian, SEXP _RatioMAD, SEXP _ExprVal)
{
  PROTECT(_RatioMedian = AS_NUMERIC(_RatioMedian));
  double *RatioMedian = NUMERIC_POINTER(_RatioMedian);
  int nGenes = INTEGER_POINTER(AS_INTEGER(GET_DIM(_RatioMedian)))[0];

  PROTECT(_RatioMAD = AS_NUMERIC(_RatioMAD));
  double *RatioMAD = NUMERIC_POINTER(_RatioMAD);
  
  PROTECT(_ExprVal = AS_NUMERIC(_ExprVal));
  double *ExprVal = NUMERIC_POINTER(_ExprVal);

  SEXP _Z;
  PROTECT(_Z = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_Z, nGenes, nGenes);
  double *Z = NUMERIC_POINTER(_Z);

  for (int i = 0; i < nGenes; i++)
    Z[i] = 0;
  
  double r;
  for (int i = 0; i < nGenes-1; i++)
  {
    for (int j = i+1; j < nGenes; j++)
    {
      if (!ISNA(ExprVal[i]) && !ISNA(ExprVal[j]) && !ISNAN(ExprVal[i]) && !ISNAN(ExprVal[j]))
      {
        r = (ExprVal[i] - ExprVal[j]) / (ExprVal[i] + ExprVal[j]);
        r = (r - RatioMedian[i+nGenes*j]) / RatioMAD[i+nGenes*j];
        
        if (!ISNAN(r))
        {
          if (r > 0) Z[i+nGenes*j] = r;
          else Z[j+nGenes*i] = -r;
        }
      }
    }
  }
  
  UNPROTECT(4);
  return(_Z);
}
