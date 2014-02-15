#include "Corbi.h"

inline int *getDepth(int *edges, int nEdges, int *core, int nGene, int *depth)
{
  for (int i = 0; i < nGene; i++)
    depth[i] = -1;

	int *out_links = (int *) R_alloc((size_t) nGene * nGene, sizeof(int));
	int *out_degree = (int *) R_alloc(nGene, sizeof(int));
	for (int i = 0; i < nGene; i++)
		out_degree[i] = 0;
  
  for (int i = 0; i < nEdges; i++)
  {
    int k = edges[i] - 1;
    out_links[(size_t) k * nGene + out_degree[k]++] = edges[i + nEdges] - 1;
  }

	bool *label_free = (bool *) R_alloc(nGene, sizeof(bool));
	int *queue = (int *) R_alloc(nGene, sizeof(int));
	int queue_head, queue_tail;

	queue_head = queue_tail = 0;

	for (int i = 0; i < nGene; i++)
	{
  	if (core[i])
    {
      queue[queue_tail++] = i;
      label_free[i] = false;
  	  depth[i] = 0;
		}
    else
    {
			label_free[i] = true;
		}
	}

  while (queue_head != queue_tail)
  {
		int k = queue[queue_head++];
		int d = depth[k];

		for (int i = 0; i < out_degree[k]; i++)
		{
			int t = out_links[(size_t) k * nGene + i];
			if (label_free[t])
			{
				queue[queue_tail++] = t;
				label_free[t] = false;
				depth[t] = d + 1;
			}
		}
	}

	return (depth);
}


SEXP NE_Depths(SEXP _Edges, SEXP _Core)
{
  PROTECT(_Edges = AS_INTEGER(_Edges));
  int *Edges = INTEGER_POINTER(_Edges);
  PROTECT(_Core = AS_LOGICAL(_Core));
  int *Core = LOGICAL_POINTER(_Core);

  SEXP _nEdges;
	PROTECT(_nEdges = GET_DIM(_Edges));
	int nEdges = INTEGER_POINTER(AS_INTEGER(_nEdges))[0];
  int nGene = length(_Core);

  SEXP _Depth;
	PROTECT(_Depth = NEW_INTEGER(nGene));
	int *Depth = INTEGER_POINTER(_Depth);

  getDepth(Edges, nEdges, Core, nGene, Depth);
  
  UNPROTECT(4);
  return (_Depth);
}

