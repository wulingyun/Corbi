#include "Corbi.h"

/* get the list element named tag, or return NULL */

SEXP getListElement(SEXP list, const char *tag)
{
	SEXP value = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), tag) == 0)
		{
			value = VECTOR_ELT(list, i);
			break;
		}
	return value;
}

/* set the list element */

void setListElement(SEXP list, int i, const char *tag, SEXP value)
{
	SEXP _names = getAttrib(list, R_NamesSymbol);
	if (_names == R_NilValue)
	{
		PROTECT(_names = NEW_STRING(length(list)));
		SET_STRING_ELT(_names, i, mkChar(tag));
		setAttrib(list, R_NamesSymbol, _names);
		UNPROTECT_PTR(_names);
	}
	else
		SET_STRING_ELT(_names, i, mkChar(tag));
	SET_VECTOR_ELT(list, i, value);
}

/* set dim of array */

void setDim2(SEXP array, int x1, int x2)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(2));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	SET_DIM(array, _dim);
	UNPROTECT_PTR(_dim);
}

void setDim3(SEXP array, int x1, int x2, int x3)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(3));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	INTEGER_POINTER(_dim)[2] = x3;
	SET_DIM(array, _dim);
	UNPROTECT_PTR(_dim);
}

void setDim4(SEXP array, int x1, int x2, int x3, int x4)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(4));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	INTEGER_POINTER(_dim)[2] = x3;
	INTEGER_POINTER(_dim)[3] = x4;
	SET_DIM(array, _dim);
	UNPROTECT_PTR(_dim);
}

int *sampleWithoutReplace(int n, int k, int *result, int *buffer)
{
  if (result == NULL)
    result = (int *) R_alloc(k, sizeof(int));
  if (buffer == NULL)
    buffer = (int *) R_alloc(n, sizeof(int));

  for (int i = 0; i < n; i++)
    buffer[i] = i;

  GetRNGstate();
  for (int i = 0; i < k; i++)
  {
    int j = n * unif_rand();
    result[i] = buffer[j];
    buffer[j] = buffer[--n];
  }
  PutRNGstate();

  return result;
}
