#include "Corbi.h"

inline int *getDepth(int *net, int *core, int n, int *depth)
{
  for (int i = 0; i < n; i++)
    depth[i] = -1;

	int *out_links = (int *) R_alloc(n * n, sizeof(int));
	int *out_degree = (int *) R_alloc(n, sizeof(int));
	for (int i = 0; i < n; i++)
	{
		out_degree[i] = 0;
		for (int j = 0; j < n; j++)
		{
			if (net[i + j * n] != 0)
			{
				out_links[i * n + out_degree[i]] = j;
				out_degree[i]++;
			}
		}
	}

	bool *label_free = (bool *) R_alloc(n, sizeof(bool));
	int *queue = (int *) R_alloc(n, sizeof(int));
	int queue_head, queue_tail;

	queue_head = queue_tail = 0;

	for (int i = 0; i < n; i++)
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
			int t = out_links[k * n + i];
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


SEXP NE_Depths(SEXP _Net, SEXP _Core)
{
  PROTECT(_Net = AS_INTEGER(_Net));
  int *Net = INTEGER_POINTER(_Net);
  PROTECT(_Core = AS_LOGICAL(_Core));
  int *Core = LOGICAL_POINTER(_Core);

  SEXP _dNet;
	PROTECT(_dNet = GET_DIM(_Net));
	int dNet = INTEGER_POINTER(AS_INTEGER(_dNet))[0];

  SEXP _Depth;
	PROTECT(_Depth = NEW_INTEGER(dNet));
	int *Depth = INTEGER_POINTER(_Depth);

  getDepth(Net, Core, dNet, Depth);
  
  UNPROTECT(4);
  return (_Depth);
}

