#include <R.h>
#include <Rdefines.h>

/* Interfaces to R */

extern "C" {
	/* Utils */
	SEXP PA_Scores(SEXP _D1, SEXP _D2, SEXP _nD1, SEXP _nD2, SEXP _afpLength);
	SEXP NQ_ShortestDistances(SEXP _W, SEXP _S);
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
