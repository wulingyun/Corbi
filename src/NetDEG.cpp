#include "Corbi.h"
#include <Rmath.h>

SEXP ND_RatioDistribution(SEXP _LogExprMatrix, SEXP _pEdge)
{
  PROTECT(_LogExprMatrix = AS_NUMERIC(_LogExprMatrix));
  double *LogExprMatrix = NUMERIC_POINTER(_LogExprMatrix);
  int *dim = INTEGER_POINTER(AS_INTEGER(GET_DIM(_LogExprMatrix)));
  int nGenes = dim[0];
  int nSamples = dim[1];

  PROTECT(_pEdge = AS_NUMERIC(_pEdge));
  double pEdge = NUMERIC_POINTER(_pEdge)[0];
  
  if (pEdge > 1) pEdge = 1;
  else if (pEdge < 0) pEdge = 0;
  double p = pEdge / 2;
  
  SEXP _LB;
  PROTECT(_LB = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_LB, nGenes, nGenes);
  double *LB = NUMERIC_POINTER(_LB);
  SetValues(_LB, LB, R_NegInf);

  double *r = (double *) R_alloc(nSamples, sizeof(double));
  double *e, m;
  int n;
  for (int i = 0; i < nGenes-1; i++)
  {
    for (int j = i+1; j < nGenes; j++)
    {
      e = LogExprMatrix;
      n = 0;
      for (int k = 0; k < nSamples; k++)
      {
        if (R_finite(e[i]) && R_finite(e[j]))
          r[n++] = e[i] - e[j];
        e += nGenes;
      }

      if (n > 0)
      {
        m = quantile(r, n, p, false);
        LB[i+nGenes*j] = m;

        m = quantile(r, n, 1-p, true);
        LB[j+nGenes*i] = -m;
      }
    }
  }

  SEXP _ratio;
  PROTECT(_ratio = NEW_LIST(2));
  SetListElement(_ratio, 0, "LB", _LB);
  SetListElement(_ratio, 1, "p.edge", _pEdge);
  
  UNPROTECT(4);
  return(_ratio);
}

SEXP ND_RatioDistribution2(SEXP _LogExprMatrix, SEXP _pEdge, SEXP _pTrim)
{
  PROTECT(_LogExprMatrix = AS_NUMERIC(_LogExprMatrix));
  double *LogExprMatrix = NUMERIC_POINTER(_LogExprMatrix);
  int *dim = INTEGER_POINTER(AS_INTEGER(GET_DIM(_LogExprMatrix)));
  int nGenes = dim[0];
  int nSamples = dim[1];
  
  PROTECT(_pEdge = AS_NUMERIC(_pEdge));
  double pEdge = NUMERIC_POINTER(_pEdge)[0];
  
  if (pEdge > 1) pEdge = 1;
  else if (pEdge < 0) pEdge = 0;

  PROTECT(_pTrim = AS_NUMERIC(_pTrim));
  double pTrim = NUMERIC_POINTER(_pTrim)[0];
  
  if (pTrim > 1) pTrim = 1;
  else if (pTrim < 0) pTrim = 0;
  
  SEXP _LB;
  PROTECT(_LB = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_LB, nGenes, nGenes);
  double *LB = NUMERIC_POINTER(_LB);
  SetValues(_LB, LB, R_NegInf);

  double *r = (double *) R_alloc(nSamples, sizeof(double));
  double *e, m, lTrim, uTrim, mTrim;
  int n, nTrim;
  for (int i = 0; i < nGenes-1; i++)
  {
    for (int j = i+1; j < nGenes; j++)
    {
      e = LogExprMatrix;
      n = 0;
      for (int k = 0; k < nSamples; k++)
      {
        if (R_finite(e[i]) && R_finite(e[j]))
          r[n++] = e[i] - e[j];
        e += nGenes;
      }
      
      if (n > 0)
      {
        lTrim = quantile(r, n, pTrim, false);
        uTrim = quantile(r, n, 1-pTrim, true);
        
        mTrim = 0;
        nTrim = 0;
        for (int k = 0; k < n; k++)
        {
          if (r[k] >= lTrim && r[k] <= uTrim)
          {
            mTrim += r[k];
            nTrim++;
          }
        }
        mTrim /= nTrim;
        
        for (int k = 0; k < n; k++)
        {
          r[k] = fabs(r[k] - mTrim);
        }
        
        m = quantile(r, n, 1-pEdge, false);
        LB[i+nGenes*j] = mTrim - m;
        LB[j+nGenes*i] = -(mTrim + m);
      }
    }
  }
  
  SEXP _ratio;
  PROTECT(_ratio = NEW_LIST(3));
  SetListElement(_ratio, 0, "LB", _LB);
  SetListElement(_ratio, 1, "p.edge", _pEdge);
  SetListElement(_ratio, 2, "p.trim", _pTrim);
  
  UNPROTECT(5);
  return(_ratio);
}

SEXP ND_DiffRatioNet(SEXP _RatioLB, SEXP _LogExprVal)
{
  PROTECT(_RatioLB = AS_NUMERIC(_RatioLB));
  double *RatioLB = NUMERIC_POINTER(_RatioLB);
  int nGenes = INTEGER_POINTER(AS_INTEGER(GET_DIM(_RatioLB)))[0];

  PROTECT(_LogExprVal = AS_NUMERIC(_LogExprVal));
  double *LogExprVal = NUMERIC_POINTER(_LogExprVal);

  int *e1 = (int *) R_alloc(nGenes * nGenes, sizeof(int));
  int *e2 = e1 + nGenes * nGenes / 2;
  int n = 0;

  double r, ei, ej, ub, lb;
  for (int i = 0; i < nGenes-1; i++)
  {
    ei = LogExprVal[i];
    if (R_finite(ei))
    {
      for (int j = i+1; j < nGenes; j++)
      {
        ej = LogExprVal[j];
        if (R_finite(ej))
        {
          r = ei - ej;
          lb = RatioLB[i+nGenes*j];
          ub = -RatioLB[j+nGenes*i];
          if (R_finite(ub) && r > ub)
          {
            e1[n] = i;
            e2[n] = j;
            n++;
          }
          else if (R_finite(lb) && r < lb)
          {
            e1[n] = j;
            e2[n] = i;
            n++;
          }
        }
      }
    }
  }

  SEXP _E1, _E2;
  PROTECT(_E1 = NEW_INTEGER(n));
  PROTECT(_E2 = NEW_INTEGER(n));
  int *E1 = INTEGER_POINTER(_E1);
  int *E2 = INTEGER_POINTER(_E2);
  
  for (int i = 0; i < n; i++)
  {
    E1[i] = e1[i] + 1;
    E2[i] = e2[i] + 1;
  }
  
  SEXP _M;
  PROTECT(_M = NEW_LIST(2));
  SetListElement(_M, 0, "i", _E1);
  SetListElement(_M, 1, "j", _E2);

  UNPROTECT(5);
  return(_M);
}
