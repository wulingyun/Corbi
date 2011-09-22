#include "Corbi.h"

SEXP NQ_ShortestDistances(SEXP _W, SEXP _S)
{
	PROTECT(_W = AS_INTEGER(_W));
	int *W = INTEGER_POINTER(_W);
	PROTECT(_S = AS_LOGICAL(_S));
	int *S = LOGICAL_POINTER(_S);

	SEXP _dim;
	PROTECT(_dim = GET_DIM(_W));
	int dim = INTEGER_POINTER(AS_INTEGER(_dim))[0];

	SEXP _D;
	PROTECT(_D = NEW_INTEGER(dim * dim));
	setDim2(_D, dim, dim);
	int *D = INTEGER_POINTER(_D);
	setValues(_D, D, -1);

	int *out_links = (int *) R_alloc(dim * dim, sizeof(int));
	int *out_degree = (int *) R_alloc(dim, sizeof(int));
	for (int i = 0; i < dim; i++)
	{
		out_degree[i] = 0;
		for (int j = 0; j < dim; j++)
		{
			if (W[i + j * dim] != 0)
			{
				out_links[i * dim + out_degree[i]] = j;
				out_degree[i]++;
			}
		}
	}

	bool *label_free = (bool *) R_alloc(dim, sizeof(bool));
	int *queue = (int *) R_alloc(dim, sizeof(int));
	int queue_head, queue_tail;

	for (int s = 0; s < dim; s++)
	{
		if (!S[s]) continue;

		for (int i = 0; i < dim; i++)
		{
			label_free[i] = true;
		}

		queue_head = queue_tail = 0;
		queue[queue_tail++] = s;
		label_free[s] = false;
		D[s + s * dim] = 0;

		while (queue_head != queue_tail)
		{
			int n = queue[queue_head++];
			int d = D[s + n * dim];

			for (int i = 0; i < out_degree[n]; i++)
			{
				int t = out_links[n * dim + i];
				if (label_free[t])
				{
					queue[queue_tail++] = t;
					label_free[t] = false;
					D[s + t * dim] = d + 1;
				}
			}
		}
	}

	UNPROTECT(4);
	return (_D);
}
