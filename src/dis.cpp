#include "Corbi.h"

SEXP dis(SEXP _W, SEXP _d1)
{
	PROTECT(_W = AS_INTEGER(_W));
	int *W = INTEGER_POINTER(_W);
	int d1 = INTEGER_POINTER(AS_INTEGER(_d1))[0];

	SEXP _D;
	PROTECT(_D = NEW_INTEGER(d1 * d1));
	setDim2(_D, d1, d1);
	int *D = INTEGER_POINTER(_D);
	setValues(_D, D, INT_MAX);

	int *out_links = (int *) R_alloc(d1 * d1, sizeof(int));
	int *out_degree = (int *) R_alloc(d1, sizeof(int));
	for (int i = 0; i < d1; i++)
	{
		out_degree[i] = 0;
		for (int j = 0; j < d1; j++)
		{
			if (W[i + j * d1] != 0)
			{
				out_links[i * d1 + out_degree[i]] = j;
				out_degree[i]++;
			}
		}
	}

	bool *label_free = (bool *) R_alloc(d1, sizeof(bool));
	int *queue = (int *) R_alloc(d1, sizeof(int));
	int queue_head, queue_tail;

	for (int s = 0; s < d1; s++)
	{
		for (int i = 0; i < d1; i++)
		{
			label_free[i] = true;
		}

		queue_head = queue_tail = 0;
		queue[queue_tail++] = s;
		label_free[s] = false;
		D[s + s * d1] = 0;

		while (queue_head != queue_tail)
		{
			int n = queue[queue_head++];
			int d = D[s + n * d1];

			for (int i = 0; i < out_degree[n]; i++)
			{
				int t = out_links[n * d1 + i];
				if (label_free[t])
				{
					queue[queue_tail++] = t;
					label_free[t] = false;
					D[s + t * d1] = d + 1;
				}
			}
		}
	}

	UNPROTECT(2);
	return (_D);
}