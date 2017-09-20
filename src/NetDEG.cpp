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