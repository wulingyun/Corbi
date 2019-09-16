#include <R.h>
#include <Rdefines.h>

/* sampling method */
/* sample without replace */

int *SampleWithoutReplace(int n, int k, int *result = NULL, int *buffer = NULL);

/* sample from discrete distribution */

int SampleFrom(int n, double *prob);

/* minimum weight spanning tree using Kruskal algorithm */

int MinSpanTree(int *tree, int nNodes, int nEdges, int *edges, double *costs, int node_index_from = 1);

/* utils for ascending ordered vector */

int Intersection(int *overlap, int *vector1, int size1, int *vector2, int size2);
int Union(int *combination, int *vector1, int size1, int *vector2, int size2);
bool Equal(int *vector1, int *vector2, int size);
void Insert(int *vector, int &size, int v);
void Remove(int *vector, int &size, int v);
