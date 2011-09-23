#include "Corbi.h"

SEXP NQ_ShortestDistances(SEXP _W, SEXP _F)
{
	PROTECT(_W = AS_NUMERIC(_W));
	double *W = NUMERIC_POINTER(_W);
	PROTECT(_F = AS_LOGICAL(_F));
	int *F = LOGICAL_POINTER(_F);

	SEXP _dW;
	PROTECT(_dW = GET_DIM(_W));
	int dW = INTEGER_POINTER(AS_INTEGER(_dW))[0];

	SEXP _D;
	PROTECT(_D = NEW_INTEGER(dW * dW));
	setDim2(_D, dW, dW);
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

	for (int s = 0; s < dW; s++)
	{
		if (!F[s]) continue;

		for (int i = 0; i < dW; i++)
		{
			label_free[i] = true;
		}

		queue_head = queue_tail = 0;
		queue[queue_tail++] = s;
		label_free[s] = false;
		D[s + s * dW] = 0;

		while (queue_head != queue_tail)
		{
			int n = queue[queue_head++];
			int d = D[s + n * dW];

			for (int i = 0; i < out_degree[n]; i++)
			{
				int t = out_links[n * dW + i];
				if (label_free[t])
				{
					queue[queue_tail++] = t;
					label_free[t] = false;
					D[s + t * dW] = d + 1;
				}
			}
		}
	}

	UNPROTECT(4);
	return (_D);
}
