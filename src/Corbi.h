#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Interfaces to R */

extern "C" {
  /* DLL Init */
  void R_init_Corbi(DllInfo *info);

	/* Utils */
  SEXP RMultiHyper(SEXP _N, SEXP _M, SEXP _K);

  SEXP PA_Scores(SEXP _D1, SEXP _D2, SEXP _nD1, SEXP _nD2, SEXP _afpLength);
	SEXP NQ_ShortestDistances(SEXP _W, SEXP _S);
  SEXP NE_Depths(SEXP _Edges, SEXP _Index, SEXP _Core);
}


/* initialize the list */
#define setValues(r, c, v) for (int i = 0; i < length(r); i++) c[i] = v

/* get/set the list element */
SEXP getListElement(SEXP list, const char *tag);
void setListElement(SEXP list, int i, const char *tag, SEXP value);

/* set dim of array */
void setDim2(SEXP array, int x1, int x2);
void setDim3(SEXP array, int x1, int x2, int x3);
void setDim4(SEXP array, int x1, int x2, int x3, int x4);

/* sampling method */
int *sampleWithoutReplace(int n, int k, int *result = NULL, int *buffer = NULL);
