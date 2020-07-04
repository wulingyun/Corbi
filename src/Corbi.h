#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include "misc.h"
#include "utils.h"
#include "quantile.h"
using namespace std;

/* Interfaces to R */

extern "C" {
  /* DLL Init */
  void R_init_Corbi(DllInfo *info);

	/* Utils */
  SEXP RMultiHyper(SEXP _N, SEXP _M, SEXP _K);
  SEXP PMultiHyper(SEXP _X, SEXP _K, SEXP _M, SEXP _W);
  SEXP PMultiNom(SEXP _X, SEXP _K, SEXP _M, SEXP _W);

	SEXP NQ_ShortestDistances(SEXP _Edges, SEXP _Index, SEXP _SourceNodes);

  SEXP PA_Scores(SEXP _D1, SEXP _D2, SEXP _nD1, SEXP _nD2, SEXP _afpLength);
  
  SEXP BS_GetSubnets(SEXP _Edges, SEXP _nNodes, SEXP _maxSize);
  SEXP BS_ExtendSubnets(SEXP _Sub1, SEXP _Sub2, SEXP _size);
  
  SEXP ND_RatioDistribution(SEXP _LogExprMatrix, SEXP _pEdge);
  SEXP ND_RatioDistributionParI(SEXP _LogExprMatrix, SEXP _pEdge, SEXP _I);
  SEXP ND_ParMerge(SEXP _SubI, SEXP _nGenes, SEXP _defV, SEXP _isSym);
  SEXP ND_RatioDistributionAB(SEXP _LogExprMatrixA, SEXP _LogExprMatrixB, SEXP _pEdge);
  SEXP ND_RatioDistributionParAiB(SEXP _LogExprMatrixAi, SEXP _LogExprMatrixB, SEXP _pEdge);
  SEXP ND_RatioDistribution2(SEXP _LogExprMatrix, SEXP _pEdge, SEXP _pTrim);
  SEXP ND_DiffRatioNet(SEXP _Dist, SEXP _LogExprVal);
  SEXP ND_RatioVariance(SEXP _LogExprMatrix);
  SEXP ND_RatioVarianceParI(SEXP _LogExprMatrix, SEXP _I);
}
