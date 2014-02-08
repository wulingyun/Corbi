#include "Corbi.h"

SEXP NE_Depths(SEXP _W, SEXP _R)
{
	PROTECT(_W = AS_NUMERIC(_W));
	double *W = NUMERIC_POINTER(_W);
	PROTECT(_R = AS_LOGICAL(_R));
	int *R = LOGICAL_POINTER(_R);

	SEXP _dW;
	PROTECT(_dW = GET_DIM(_W));
	int dW = INTEGER_POINTER(AS_INTEGER(_dW))[0];

	SEXP _D;
	PROTECT(_D = NEW_INTEGER(dW));
	int *D = INTEGER_POINTER(_D);
	setValues(_D, D, -1);

	int *out_links = (int *) R_alloc(dW * dW, sizeof(int));
	int *out_degree = (int *) R_alloc(dW, sizeof(int));
	for (int i = 0; i < dW; i++)
	{
		out_degree[i] = 0;
		for (int j = 0; j < dW; j++)
		{
			if (W[i + j * dW] != 0)
			{
				out_links[i * dW + out_degree[i]] = j;
				out_degree[i]++;
			}
		}
	}

	bool *label_free = (bool *) R_alloc(dW, sizeof(bool));
	int *queue = (int *) R_alloc(dW, sizeof(int));
	int queue_head, queue_tail;

	queue_head = queue_tail = 0;

	for (int i = 0; i < dW; i++)
	{
  	if (R[i])
    {
      queue[queue_tail++] = i;
      label_free[i] = false;
  	  D[i] = 0;
		}
    else
    {
			label_free[i] = true;
		}
	}

  while (queue_head != queue_tail)
  {
		int n = queue[queue_head++];
		int d = D[n];

		for (int i = 0; i < out_degree[n]; i++)
		{
			int t = out_links[n * dW + i];
			if (label_free[t])
			{
				queue[queue_tail++] = t;
				label_free[t] = false;
				D[t] = d + 1;
			}
		}
	}

	UNPROTECT(4);
	return (_D);
}
