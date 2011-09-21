#include "Corbi.h"
#include <Rmath.h>

inline double _afp_score(double *D1, double *D2, int nD1, int nD2, int afpLength, int i, int j)
{
	double d = 0;
	for (int k1 = 0; k1 < afpLength; k1++)
	{
		for (int k2 = k1+1; k2 < afpLength; k2++)
		{
			d += R_pow_di(D1[i+k1 + nD1 * (i+k2)] - D2[j+k1 + nD2 * (j+k2)], 2);
		}
	}
	d /= afpLength * (afpLength - 1) / 2;
	return (d);
}

inline double _afp_dist(double *D1, double *D2, int nD1, int nD2, int afpLength, int i1, int i2, int j1, int j2)
{
	double d = 0;
	for (int k = 0; k < afpLength; k++)
	{
		d += R_pow_di(D1[i1+k + nD1 * (i2+k)] - D2[j1+k + nD2 * (j2+k)], 2);
	}
	d /= afpLength;
	return (d);
}

template <class T>
inline const T& min(const T& a, const T& b)
{
	return a < b ? a : b;
}

template <class T>
inline const T& max(const T& a, const T& b)
{
	return a > b ? a : b;
}

SEXP PA_Scores(SEXP _D1, SEXP _D2, SEXP _nD1, SEXP _nD2, SEXP _afpLength)
{
	PROTECT(_D1 = AS_NUMERIC(_D1));
	PROTECT(_D2 = AS_NUMERIC(_D2));
	double *D1 = NUMERIC_POINTER(_D1);
	double *D2 = NUMERIC_POINTER(_D2);

	int nD1 = INTEGER_POINTER(AS_INTEGER(_nD1))[0];
	int nD2 = INTEGER_POINTER(AS_INTEGER(_nD2))[0];
	int afpLength = INTEGER_POINTER(AS_INTEGER(_afpLength))[0];
	int nAFP1 = nD1 - afpLength + 1;
	int nAFP2 = nD2 - afpLength + 1;

	/* Node score */

	SEXP _nodeScore;
	PROTECT(_nodeScore = NEW_NUMERIC(nD1 * nD2));
	setDim2(_nodeScore, nD1, nD2);
	double *nodeScore = NUMERIC_POINTER(_nodeScore);

	double *afpScore = (double *) R_alloc(nAFP1 * nAFP2, sizeof(double));
	for (int i = 0; i < nAFP1; i++)
	{
		for (int j = 0; j < nAFP2; j++)
		{
			afpScore[i + nAFP1 * j] = _afp_score(D1, D2, nD1, nD2, afpLength, i, j);
		}
	}

	for (int i = 0; i < nD1; i++)
	{
		int kMin0 = max(0, i - nAFP1 + 1);
		int kMax0 = min(afpLength, i + 1);
		for (int j = 0; j < nD2; j++)
		{
			int kMin = max(kMin0, j - nAFP2 + 1);
			int kMax = min(kMax0, j + 1);
			double d = 1e100;
			for (int k = kMin; k < kMax; k++)
				d = min(d, afpScore[i-k + nAFP1 * (j-k)]);
			nodeScore[i + nD1 * j] = d;
		}
	}
	
	/* Edge score */

	SEXP _edgeScore;
	PROTECT(_edgeScore = NEW_NUMERIC(nD2 * nD2 * (nD1 - 1)));
	setDim3(_edgeScore, nD2, nD2, (nD1 - 1));
	double *edgeScore = NUMERIC_POINTER(_edgeScore);

	double *afpDist = (double *) R_alloc(nAFP2 * nAFP2 * afpLength * (nAFP1 - 1), sizeof(double));
	for (int i1 = 0; i1 < (nAFP1 - 1); i1++)
	{
		int kMax = min(afpLength, nAFP1 - 1 - i1);
		for (int k = 0; k < kMax; k++)
		{
			int i2 = i1 + 1 + k;
			for (int j1 = 0; j1 < nAFP2; j1++)
			{
				for (int j2 = 0; j2 < nAFP2; j2++)
				{
					afpDist[j1 + nAFP2 * (j2 + nAFP2 * (k + afpLength * i1))] = _afp_dist(D1, D2, nD1, nD2, afpLength, i1, i2, j1, j2);
				}
			}
		}
	}

	for (int i1 = 0; i1 < (nD1 - 1); i1++)
	{
		int k1Min0 = max(0, i1 - nAFP1 + 1);
		int k1Max0 = min(afpLength, i1 + 1);

		int i2 = i1 + 1;
		int k2Min0 = max(0, i2 - nAFP1 + 1);
		int k2Max0 = min(afpLength, i2 + 1);
		for (int j1 = 0; j1 < nD2; j1++)
		{
			int k1Min = max(k1Min0, j1 - nAFP2 + 1);
			int k1Max = min(k1Max0, j1 + 1);
			for (int j2 = 0; j2 < nD2; j2++)
			{
				int k2Min = max(k2Min0, j2 - nAFP2 + 1);
				int k2Max = min(k2Max0, j2 + 1);

				double d = 1e100;
				if (j1 < j2 - 1 && k2Min == 0)
				{
					for (int k = k1Min; k < k1Max; k++)
					{
						d = min(d, afpDist[j1-k + nAFP2 * (j2 + nAFP2 * (k + afpLength * (i1 - k)))]);
					}
				}
				else if (j1 == j2 - 1)
				{
					for (int k = k1Min; k < k2Max - 1; k++)
					{
						d = min(d, afpScore[i1-k + nAFP1 * (j1-k)]);
					}
				}
				edgeScore[j1 + nD2 * (j2 + nD2 * i1)] = d;
			}
		}
	}

	SEXP _scores;
	PROTECT(_scores = NEW_LIST(2));
	setListElement(_scores, 0, "node.score", _nodeScore);
	setListElement(_scores, 1, "edge.score", _edgeScore);

	UNPROTECT(5);
	return (_scores);
}
