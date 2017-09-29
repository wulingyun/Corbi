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

SEXP ND_RatioDistribution(SEXP _ExprVal, SEXP _pEdge)
{
  PROTECT(_ExprVal = AS_NUMERIC(_ExprVal));
  double *ExprVal = NUMERIC_POINTER(_ExprVal);
  int *dim = INTEGER_POINTER(AS_INTEGER(GET_DIM(_ExprVal)));
  int nGenes = dim[0];
  int nSamples = dim[1];

  PROTECT(_pEdge = AS_NUMERIC(_pEdge));
  double pEdge = NUMERIC_POINTER(_pEdge)[0];
  
  if (pEdge > 1) pEdge = 1;
  if (pEdge < 0) pEdge = 0;
  double p = pEdge / 2;
  
  SEXP _LB;
  PROTECT(_LB = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_LB, nGenes, nGenes);
  double *LB = NUMERIC_POINTER(_LB);
  
  SEXP _UB;
  PROTECT(_UB = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_UB, nGenes, nGenes);
  double *UB = NUMERIC_POINTER(_UB);

  for (int i = 0; i < nGenes * nGenes; i++)
  {
    UB[i] = 0;
    LB[i] = 0;
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
          r[n++] = log(e[i]) - log(e[j]);
          // r[n++] = (e[i] - e[j]) / (e[i] + e[j]);
        e += nGenes;
      }

      m = quantile(r, n, p, false);
      LB[i+nGenes*j] = m;
      UB[j+nGenes*i] = -m;

      m = quantile(r, n, 1-p, true);
      UB[i+nGenes*j] = m;
      LB[j+nGenes*i] = -m;
    }
  }

  SEXP _ratio;
  PROTECT(_ratio = NEW_LIST(3));
  SetListElement(_ratio, 0, "LB", _LB);
  SetListElement(_ratio, 1, "UB", _UB);
  SetListElement(_ratio, 2, "p.edge", _pEdge);
  
  UNPROTECT(5);
  return(_ratio);
}

SEXP ND_DiffRatioNet(SEXP _RatioLB, SEXP _RatioUB, SEXP _ExprVal)
{
  PROTECT(_RatioLB = AS_NUMERIC(_RatioLB));
  double *RatioLB = NUMERIC_POINTER(_RatioLB);
  int nGenes = INTEGER_POINTER(AS_INTEGER(GET_DIM(_RatioLB)))[0];

  PROTECT(_RatioUB = AS_NUMERIC(_RatioUB));
  double *RatioUB = NUMERIC_POINTER(_RatioUB);
  
  PROTECT(_ExprVal = AS_NUMERIC(_ExprVal));
  double *ExprVal = NUMERIC_POINTER(_ExprVal);

  SEXP _M;
  PROTECT(_M = NEW_INTEGER(nGenes * nGenes));
  SetDim2(_M, nGenes, nGenes);
  int *M = INTEGER_POINTER(_M);

  for (int i = 0; i < nGenes*nGenes; i++)
    M[i] = 0;
  
  double r;
  for (int i = 0; i < nGenes-1; i++)
  {
    for (int j = i+1; j < nGenes; j++)
    {
      if (!ISNA(ExprVal[i]) && !ISNA(ExprVal[j]) && !ISNAN(ExprVal[i]) && !ISNAN(ExprVal[j]))
      {
        r = log(ExprVal[i]) - log(ExprVal[j]);
        // r = (ExprVal[i] - ExprVal[j]) / (ExprVal[i] + ExprVal[j]);
        
        if (!ISNAN(r))
        {
          if (r > RatioUB[i+nGenes*j]) M[i+nGenes*j] = 1;
          else if (r < RatioLB[i+nGenes*j]) M[j+nGenes*i] = 1;
        }
      }
    }
  }
  
  UNPROTECT(4);
  return(_M);
}
