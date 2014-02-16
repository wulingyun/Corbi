#include "Corbi.h"
#include <Rmath.h>

SEXP RMultiHyper(SEXP _N, SEXP _M, SEXP _K)
{
  int N = INTEGER_POINTER(AS_INTEGER(_N))[0];
  int K = INTEGER_POINTER(AS_INTEGER(_K))[0];

  PROTECT(_M = AS_INTEGER(_M));
  int *M = INTEGER_POINTER(_M);
  int nM = length(_M);

  SEXP _X;
  PROTECT(_X = NEW_INTEGER(nM * N));
  setDim2(_X, nM, N);
  int *X = INTEGER_POINTER(_X);

  int sumM = 0;
  for (int i = 0; i < nM; i++)
      sumM += M[i];

  GetRNGstate();
  for (int n = 0; n < N; n++)
  {
    int restM = sumM;
    int restK = K;
    for (int i = 1; i < nM; i++)
    {
      restM -= M[i];
      X[i] = rhyper(M[i], restM, restK);
      restK -= X[i];
    }
    X[0] = restK;
    X += nM;
  }
  PutRNGstate();
  
  UNPROTECT(2);
  return (_X);
}
