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
      e = LogExprMatrix;
      n = 0;
      for (int k = 0; k < nSamples; k++)
      {
        if (!ISNA(e[i]) && !ISNA(e[j]) && !ISNAN(e[i]) && !ISNAN(e[j]))
          r[n++] = e[i] - e[j];
        e += nGenes;
      }

      if (n > 0)
      {
        m = quantile(r, n, p, false);
        LB[i+nGenes*j] = m;
        UB[j+nGenes*i] = -m;
        
        m = quantile(r, n, 1-p, true);
        UB[i+nGenes*j] = m;
        LB[j+nGenes*i] = -m;
      }
      else
      {
        LB[i+nGenes*j] = R_NegInf;
        LB[j+nGenes*i] = R_NegInf;
        UB[j+nGenes*i] = R_PosInf;
        UB[i+nGenes*j] = R_PosInf;
      }
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

SEXP ND_DiffRatioNet(SEXP _RatioLB, SEXP _RatioUB, SEXP _LogExprVal)
{
  PROTECT(_RatioLB = AS_NUMERIC(_RatioLB));
  double *RatioLB = NUMERIC_POINTER(_RatioLB);
  int nGenes = INTEGER_POINTER(AS_INTEGER(GET_DIM(_RatioLB)))[0];

  PROTECT(_RatioUB = AS_NUMERIC(_RatioUB));
  double *RatioUB = NUMERIC_POINTER(_RatioUB);
  
  PROTECT(_LogExprVal = AS_NUMERIC(_LogExprVal));
  double *LogExprVal = NUMERIC_POINTER(_LogExprVal);

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
      double ei, ej, ub, lb;
      ei = LogExprVal[i];
      ej = LogExprVal[j];
      ub = RatioUB[i+nGenes*j];
      lb = RatioLB[i+nGenes*j];
      if (!ISNA(ei) && !ISNA(ej) && !ISNAN(ei) && !ISNAN(ej))
      {
        r = ei - ej;
        if (!ISNA(ub) && !ISNAN(ub) && r > ub) M[i+nGenes*j] = 1;
        else if (!ISNA(lb) && !ISNAN(lb) && r < lb) M[j+nGenes*i] = 1;
      }
    }
  }
  
  UNPROTECT(4);
  return(_M);
}
