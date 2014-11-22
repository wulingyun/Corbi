#include "Corbi.h"

/* get the list element named tag, or return NULL */

SEXP GetListElement(SEXP list, const char *tag)
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

void SetListElement(SEXP list, int i, const char *tag, SEXP value)
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

void SetDim2(SEXP array, int x1, int x2)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(2));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	SET_DIM(array, _dim);
	UNPROTECT_PTR(_dim);
}

void SetDim3(SEXP array, int x1, int x2, int x3)
{
	SEXP _dim;
	PROTECT(_dim = NEW_INTEGER(3));
	INTEGER_POINTER(_dim)[0] = x1;
	INTEGER_POINTER(_dim)[1] = x2;
	INTEGER_POINTER(_dim)[2] = x3;
	SET_DIM(array, _dim);
	UNPROTECT_PTR(_dim);
}

void SetDim4(SEXP array, int x1, int x2, int x3, int x4)
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

/* sampling method */

/* sample without replace */

int *SampleWithoutReplace(int n, int k, int *result, int *buffer)
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

/* sample from discret distribution */

int SampleFrom(int n, double *prob)
{
  int select = n-1;
	double cutoff = unif_rand();
	double cumulativeProb = 0;
	for (int i = 0; i < n; i++)
	{
		cumulativeProb += prob[i];
		if (cumulativeProb > cutoff)
		{
			select = i;
			break;
		}
	}
	return select;
}

/* minimum weight spanning tree using Kruskal algorithm */

int MinSpanTree(int *tree, int nNodes, int nEdges, int *edges, double *costs, int node_index_from)
{
	int *index = (int *) C_allocVector<int>(nEdges);
	for (int i = 0; i < nEdges; i++)
	{
		tree[i] = 0;
		index[i] = i;
	}
	rsort_with_index(costs, index, nEdges);

	int *label = (int *) C_allocVector<int>(nNodes);
	for (int i = 0; i < nNodes; i++)
		label[i] = i;

	int n = 0, n1, n2;
	for (int i = 0; i < nEdges; i++)
	{
		n1 = label[edges[index[i]] - node_index_from];
		n2 = label[edges[index[i] + nEdges] - node_index_from];
		if (n1 != n2)
		{
			for (int j = 0; j < nNodes; j++)
				if (label[j] == n2)
					label[j] = n1;
			tree[index[i]] = 1;
			if (++n >= nNodes - 1)
				break;
		}
	}

	C_freeVector(index);
	C_freeVector(label);

	return n;
}

/* Get intersection of two ascending ordered vectors */

int Intersection(int *overlap, int *vector1, int size1, int *vector2, int size2)
{
	int n, i1, i2;
	if (vector1[0] > vector2[size2-1] || vector2[0] > vector1[size1-1])
		return 0;
	n = i1 = i2 = 0;
	while (i1 < size1 && i2 < size2)
	{
		if (vector1[i1] == vector2[i2])
		{
			overlap[n++] = vector1[i1++];
			i2++;
		}
		else if (vector1[i1] < vector2[i2])
			i1++;
		else
			i2++;
	}
	return n;
}

/* Get union of two ascending ordered vectors */

int Union(int *combination, int *vector1, int size1, int *vector2, int size2)
{
  int n, i1, i2;
	if (vector1[0] > vector2[size2-1])
  {
    for (int i = 0; i < size2; i++)
      combination[i] = vector2[i];
    for (int i = 0; i < size1; i++)
      combination[size2+i] = vector1[i];
    n = size1 + size2;
  }
  else if (vector2[0] > vector1[size1-1])
  {
    for (int i = 0; i < size1; i++)
      combination[i] = vector1[i];
    for (int i = 0; i < size2; i++)
      combination[size1+i] = vector2[i];
    n = size1 + size2;
  }
	else
  {
    n = i1 = i2 = 0;
  	while (i1 < size1 && i2 < size2)
  	{
  		if (vector1[i1] == vector2[i2])
  		{
  			combination[n++] = vector1[i1++];
  			i2++;
  		}
  		else if (vector1[i1] < vector2[i2])
        combination[n++] = vector1[i1++];
  		else
        combination[n++] = vector2[i2++];
  	}
    while (i1 < size1)
      combination[n++] = vector1[i1++];
    while (i2 < size2)
      combination[n++] = vector2[i2++];
  }
	return n;
}

/* Check equal of two ascending ordered vectors */

bool Equal(int *vector1, int *vector2, int size)
{
  bool equal = true;
  for (int i = 0; i < size; i++)
  {
    if (vector1[i] != vector2[i])
    {
      equal = false;
      break;
    }
  }
  return equal;
}

/* Insert element to ascending ordered vector */

void Insert(int *vector, int &size, int v)
{
	int k = size;
	for (int i = 0; i < size; i++)
	{
		if (vector[i] > v)
		{
			for (int j = size; j > i; j--)
				vector[j] = vector[j-1];
			k = i;
			break;
		}
	}
	vector[k] = v;
	size++;
}

/* Remove element from ascending ordered vector */

void Remove(int *vector, int &size, int v)
{

	for (int i = 0; i < size; i++)
	{
		if (vector[i] == v)
		{
			for (int j = i; j < size-1; j++)
				vector[j] = vector[j+1];
			size--;
			break;
		}
	}
}
