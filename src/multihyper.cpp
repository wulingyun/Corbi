#include "Corbi.h"
#include <Rmath.h>

SEXP RMultiHyper(SEXP _N, SEXP _M, SEXP _K)
{
  int N = INTEGER_POINTER(AS_INTEGER(_N))[0];
  int K = INTEGER_POINTER(AS_INTEGER(_K))[0];

  PROTECT(_M = AS_INTEGER(_M));
  int *M = INTEGER_POINTER(_M);
  int len = length(_M);

  SEXP _X;
  PROTECT(_X = NEW_INTEGER(len * N));
  setDim2(_X, len, N);
  int *X = INTEGER_POINTER(_X);

  GetRNGstate();
  for (int n = 0; n < N; n++)
  {
    int nM = 0;
    for (int i = 0; i < len; i++)
      nM += M[i];
    int nK = K;
    for (int i = 1; i < len; i++)
    {
      nM -= M[i];
      X[i] = rhyper(M[i], nM, nK);
      nK -= X[i];
    }
    X[0] = nK;
    X += len;
  }
  PutRNGstate();
  
  UNPROTECT(2);
  return (_X);
}
