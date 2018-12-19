#include "Corbi.h"

class Tree {
public:
  int value;
  Tree *child;
  Tree *next;
  
  Tree(int v = -1);
  ~Tree();
  
  bool Add(int *vector, int size);
  bool Find(int *vector, int size);
  int *Export(int *matrix, int nRow, int size);
};

Tree::Tree(int v)
{
  value = v;
  child = NULL;
  next = NULL;
}

Tree::~Tree()
{
  delete child;
  delete next;
}

bool Tree::Add(int *vector, int size)
{
  bool added = false;
  Tree *p = this, *p0, *p1;
  for (int i = 0; i < size; i++)
  {
    p0 = p;
    while (p->value < vector[i] && p->next != NULL)
    {
      p0 = p;
      p = p->next;
    }
    if (p->value < vector[i])
    {
      added = true;
      p1 = new Tree(vector[i]);
      p->next = p1;
    }
    else if (p->value > vector[i])
    {
      added = true;
      p1 = new Tree(vector[i]);
      p0->next = p1;
      p1->next = p;
    }
    else
      p1 = p;
    if (i < (size-1) && p1->child == NULL)
    {
      p = new Tree();
      p1->child = p;
    }
    else
      p = p1->child;
  }
  return added;
}

bool Tree::Find(int *vector, int size)
{
  bool found = false;
  Tree *p = this;
  for (int i = 0; i < size; i++)
  {
    while (p->value < vector[i] && p->next != NULL)
      p = p->next;
    if (p->value == vector[i])
    {
      if (p->child != NULL)
        p = p->child;
      else if (i == (size-1))
        found = true;
      else
        break;
    }
    else
      break;
  }
  return found;
}

int *Tree::Export(int *matrix, int nRow, int size)
{
  Tree **p = (Tree **) Calloc(size, Tree *);
  
  int n = 0, i = 0;
  p[0] = this;
  while (i >= 0)
  {
    while (p[i]->child != NULL)
    {
      p[i+1] = p[i]->child;
      i++;
    }
    if (p[i]->value >= 0)
    {
      for (int k = 0; k < size; k++)
        matrix[k * nRow + n] = p[k]->value;
      n++;
    }
    while (i >= 0 && p[i]->next == NULL)
    {
      i--;
    }
    if (i >= 0)
      p[i] = p[i]->next;
  }

  Free(p);
  return matrix;
}

SEXP extend(int *Sub1, int *Sub2, int n1, int n2, int s1, int s2, int size)
{
  int **sub1 = C_allocArray<int>(n1, s1);
  for (int i = 0; i < n1; i++)
    for (int k = 0; k < s1; k++)
      sub1[i][k] = Sub1[k * n1 + i];

  int **sub2;
  if (Sub2 != Sub1)
  {
    sub2 = C_allocArray<int>(n2, s2);
    for (int i = 0; i < n2; i++)
      for (int k = 0; k < s2; k++)
        sub2[i][k] = Sub2[k * n2 + i];
  }
  else
  {
    sub2 = sub1;
  }

  Tree tree;
  int *temp = C_allocVector<int>(s1 + s2);

  int nSub = 0, j0;
  for (int i = 0; i < n1; i++)
  {
    if (Sub2 != Sub1)
      j0 = 0;
    else
      j0 = i+1;
    for (int j = j0; j < n2; j++)
    {
      if (Union(temp, sub1[i], s1, sub2[j], s2) == size && tree.Add(temp, size))
          nSub++;
    }
  }

  SEXP _sub;
  PROTECT(_sub = NEW_INTEGER(nSub * size));
  SetDim2(_sub, nSub, size);
  int *sub = INTEGER_POINTER(_sub);
  tree.Export(sub, nSub, size);

  C_freeArray<int, 2>(sub1);
  if (Sub2 != Sub1) C_freeArray<int, 2>(sub2);
  C_freeVector(temp);

  UNPROTECT(1);
  return(_sub);
}

SEXP BS_GetSubnets(SEXP _Edges, SEXP _nNodes, SEXP _maxSize)
{
  PROTECT(_Edges = AS_INTEGER(_Edges));
  int *Edges = INTEGER_POINTER(_Edges);
  int nEdges = INTEGER_POINTER(GET_DIM(_Edges))[0];

  int nNodes = INTEGER_POINTER(AS_INTEGER(_nNodes))[0];
	int maxSize = INTEGER_POINTER(AS_INTEGER(_maxSize))[0];
  
  if (maxSize < 2) maxSize = 2;
  if (maxSize > nNodes) maxSize = nNodes;

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
    PROTECT(_temp = extend(Subnets[i-1], Subnets[1], n, nEdges, i, 2, i+1));
  	Subnets[i] = INTEGER_POINTER(_temp);
  	SET_VECTOR_ELT(_Subnets, i, _temp);
    n = INTEGER_POINTER(GET_DIM(_temp))[0];
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
