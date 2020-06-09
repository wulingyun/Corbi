#include "Corbi.h"

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

SEXP ND_RatioDistributionParI(SEXP _LogExprMatrix, SEXP _pEdge, SEXP _I)
{
  PROTECT(_LogExprMatrix = AS_NUMERIC(_LogExprMatrix));
  double *LogExprMatrix = NUMERIC_POINTER(_LogExprMatrix);
  int *dim = INTEGER_POINTER(AS_INTEGER(GET_DIM(_LogExprMatrix)));
  int nGenes = dim[0];
  int nSamples = dim[1];

  PROTECT(_pEdge = AS_NUMERIC(_pEdge));
  double pEdge = NUMERIC_POINTER(_pEdge)[0];

  PROTECT(_I = AS_INTEGER(_I));
  int I = INTEGER_POINTER(_I)[0];
    
  if (pEdge > 1) pEdge = 1;
  else if (pEdge < 0) pEdge = 0;
  double p = pEdge / 2;
  
  int nLB = nGenes - I;
  SEXP _LB;
  PROTECT(_LB = NEW_NUMERIC(nGenes * 2));
  SetDim2(_LB, nGenes, 2);
  double *LB = NUMERIC_POINTER(_LB);
  SetValues(_LB, LB, R_NegInf);

  double *r = (double *) R_alloc(nSamples, sizeof(double));
  double *e, m;
  int n;

  int II[2], KK;
  II[0] = I - 1;
  II[1] = nGenes - I - 1;
  if (II[0] < II[1]) KK = 2;
  else KK = 1;
  
  for (int ii = 0; ii < KK; ii++)
  {
    int i = II[ii];
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
        LB[j-i-1 + nLB*ii] = m;
  
        m = quantile(r, n, 1-p, true);
        LB[j-i-1 + nLB*ii + nGenes] = -m;
      }
    }
  }

  UNPROTECT(4);
  return(_LB);
}

SEXP ND_RatioDistributionParM(SEXP _DistI, SEXP _nGenes)
{
  PROTECT(_nGenes = AS_INTEGER(_nGenes));
  int nGenes = INTEGER_POINTER(_nGenes)[0];

  SEXP _LB;
  PROTECT(_LB = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_LB, nGenes, nGenes);
  double *LB = NUMERIC_POINTER(_LB);
  SetValues(_LB, LB, R_NegInf);

  for (int I = 0; I < (nGenes/2); I++)
  {
    SEXP _M;
    PROTECT(_M = AS_NUMERIC(GetListElement(_DistI, I)));
    double *m = NUMERIC_POINTER(_M);
    
    int II[2], KK, nLB;
    II[0] = I;
    II[1] = nGenes - I - 2;
    if (II[0] < II[1]) KK = 2;
    else KK = 1;
    nLB = nGenes - I - 1;
  
    for (int ii = 0; ii < KK; ii++)
    {
      int i = II[ii];
      for (int j = i+1; j < nGenes; j++)
      {
        LB[i+nGenes*j] = m[j-i-1 + nLB*ii];
        LB[j+nGenes*i] = m[j-i-1 + nLB*ii + nGenes];
      }
    }
    UNPROTECT(1);
  }

  SEXP _ratio;
  PROTECT(_ratio = NEW_LIST(1));
  SetListElement(_ratio, 0, "LB", _LB);

  UNPROTECT(3);
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

SEXP ND_RatioVariance(SEXP _LogExprMatrix)
{
  PROTECT(_LogExprMatrix = AS_NUMERIC(_LogExprMatrix));
  double *LogExprMatrix = NUMERIC_POINTER(_LogExprMatrix);
  int *dim = INTEGER_POINTER(AS_INTEGER(GET_DIM(_LogExprMatrix)));
  int nGenes = dim[0];
  int nSamples = dim[1];

  SEXP _Var;
  PROTECT(_Var = NEW_NUMERIC(nGenes * nGenes));
  SetDim2(_Var, nGenes, nGenes);
  double *Var = NUMERIC_POINTER(_Var);
  SetValues(_Var, Var, 0.0);
  
  double *r = (double *) R_alloc(nSamples, sizeof(double));
  double *e, m, v;
  int n;
  for (int i = 0; i < nGenes-1; i++)
  {
    for (int j = i+1; j < nGenes; j++)
    {
      e = LogExprMatrix;
      n = 0;
      m = v = 0;
      for (int k = 0; k < nSamples; k++)
      {
        if (R_finite(e[i]) && R_finite(e[j]))
        {
          r[n] = e[i] - e[j];
          m += r[n];
          v += r[n] * r[n];
          n++;
        }
        e += nGenes;
      }
      
      if (n > 0)
      {
        v = (v - m * m / n)  / (n - 1);
        Var[i+nGenes*j] = v;
        Var[j+nGenes*i] = v;
      }
    }
  }

  UNPROTECT(2);
  return(_Var);
}
