#include "Corbi.h"

SEXP extend(int *Sub1, int *Sub2, int n1, int n2, int s1, int s2, int size)
{
  int **sub1 = C_allocArray<int>(n1, s1);
  for (int i = 0; i < n1; i++)
    for (int k = 0; k < s1; k++)
      sub1[i][k] = Sub1[k * n1 + i];

  int **sub2 = C_allocArray<int>(n2, s2);
  for (int i = 0; i < n2; i++)
    for (int k = 0; k < s2; k++)
      sub2[i][k] = Sub2[k * n2 + i];

  int nTemp = n1 * n2;
  int **temp = C_allocArray<int>(nTemp, size);
  int *v1 = C_allocVector<int>(n1 + n2);
  int *v2 = C_allocVector<int>(n1 + n2);

  int nSub = 0, t1, t2;
  for (int i = 0; i < n1; i++)
  {
    for (int j = 0; j < n2; j++)
    {
      t1 = Union(v1, sub1[i], s1, sub2[j], s2);
      if (t1 == size)
      {
        bool duplicated = false;
        for (int k = 0; k < nSub; k++)
        {
          t2 = Intersection(v2, v1, t1, temp[k], size);
          if (t2 == size)
          {
            duplicated = true;
            break;
          }
        }
        if (!duplicated)
        {
          for (int k = 0; k < size; k++)
            temp[nSub][k] = v1[k];
          nSub++;
        }
      }
    }
  }

  SEXP _sub;
  PROTECT(_sub = NEW_INTEGER(nSub * size));
  SetDim2(_sub, nSub, size);
  int *sub = INTEGER_POINTER(_sub);
  for (int i = 0; i < nSub; i++)
    for (int k = 0; k < size; k++)
      sub[k * nSub + i] = temp[i][k];

  C_freeArray<int, 2>(sub1);
  C_freeArray<int, 2>(sub2);
  C_freeArray<int, 2>(temp);
  C_freeVector(v1);
  C_freeVector(v2);

  UNPROTECT(1);
  return(_sub);
}

SEXP BS_GetSubnets(SEXP _Edges, SEXP _nNodes, SEXP _maxSize)
{
  PROTECT(_Edges = AS_INTEGER(_Edges));
  int *Edges = INTEGER_POINTER(_Edges);
  int nEdges = INTEGER_POINTER(AS_INTEGER(GET_DIM(_Edges)))[0];

  int nNodes = INTEGER_POINTER(AS_INTEGER(_nNodes))[0];
	int maxSize = INTEGER_POINTER(AS_INTEGER(_maxSize))[0];
  
  if (maxSize < 2) maxSize = 2;

	SEXP _Subnets;
	PROTECT(_Subnets = NEW_LIST(maxSize));
  int **Subnets = (int **) R_alloc(maxSize, sizeof(int *));
  SEXP _temp;

  // subnets of size 1
	PROTECT(_temp = NEW_INTEGER(nNodes));
	SetDim2(_temp, nNodes, 1);
  SET_VECTOR_ELT(_Subnets, 0, _temp);
  Subnets[0] = INTEGER_POINTER(_temp);
  for (int i = 0; i < nNodes; i++)
  {
    Subnets[0][i] = i+1;
  }

  // subnets of size 2
  PROTECT(_temp = NEW_INTEGER(nEdges * 2));
	SetDim2(_temp, nEdges, 2);
	SET_VECTOR_ELT(_Subnets, 1, _temp);
  Subnets[1] = INTEGER_POINTER(_temp);
  for (int i = 0; i < nEdges*2; i++)
  {
    Subnets[1][i] = Edges[i];
  }

  // subnets of size > 2
  int n = nEdges;
  for (int i = 2; i < maxSize; i++)
  {
    PROTECT(_temp = extend(Subnets[i-1], Subnets[i-1], n, n, i, i, i+1));
  	Subnets[i] = INTEGER_POINTER(_temp);
  	SET_VECTOR_ELT(_Subnets, i, _temp);
    n = INTEGER_POINTER(AS_INTEGER(GET_DIM(_temp)))[0];
  }

  UNPROTECT(maxSize + 2);
	return (_Subnets);
}

SEXP BS_ExtendSubnets(SEXP _Sub1, SEXP _Sub2, SEXP _size)
{
  PROTECT(_Sub1 = AS_INTEGER(_Sub1));
  int *Sub1 = INTEGER_POINTER(_Sub1);
  PROTECT(_Sub2 = AS_INTEGER(_Sub2));
  int *Sub2 = INTEGER_POINTER(_Sub2);
  
  SEXP _nSub1, _nSub2;
  PROTECT(_nSub1 = GET_DIM(_Sub1));
  PROTECT(_nSub2 = GET_DIM(_Sub2));
  int n1 = INTEGER_POINTER(AS_INTEGER(_nSub1))[0];
  int s1 = INTEGER_POINTER(AS_INTEGER(_nSub1))[1];
  int n2 = INTEGER_POINTER(AS_INTEGER(_nSub2))[0];
  int s2 = INTEGER_POINTER(AS_INTEGER(_nSub2))[1];

  int size = INTEGER_POINTER(AS_INTEGER(_size))[0];

  SEXP _Subnets = extend(Sub1, Sub2, n1, n2, s1, s2, size);

  UNPROTECT(4);
  return(_Subnets);
}
