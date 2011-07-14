#include "Corbi.h"
#include <Rmath.h>

inline double _afp_score(double *D1, double *D2, int nD1, int nD2, int afpLength, int i, int j)
{
	double d = 0;
	for (int k1 = 0; k1 < afpLength; k1++)
		for (int k2 = k1+1; k2 < afpLength; k2++)
		{
			d += R_pow_di(D1[i+k1-1 + nD1 * (i+k2-1)] - D2[j+k1-1 + nD2 * (j+k2-1)], 2);
		}
	d = sqrt(d/(afpLength * (afpLength - 1) / 2));
	return (d);
}

SEXP AFP_Score(SEXP _D1, SEXP _D2, SEXP _nD1, SEXP _nD2, SEXP _afpLength)
{
	PROTECT(_D1 = AS_NUMERIC(_D1));
	PROTECT(_D2 = AS_NUMERIC(_D2));
	double *D1 = NUMERIC_POINTER(_D1);
	double *D2 = NUMERIC_POINTER(_D2);

	int nD1 = INTEGER_POINTER(AS_INTEGER(_nD1))[0];
	int nD2 = INTEGER_POINTER(AS_INTEGER(_nD2))[0];
	int afpLength = INTEGER_POINTER(AS_INTEGER(_afpLength))[0];

	int n1 = nD1-afpLength+1, n2 = nD2-afpLength+1;

	SEXP _mScore;
	PROTECT(_mScore = NEW_NUMERIC(n1 * n2));
	setDim2(_mScore, n1, n2);
	double *mScore = NUMERIC_POINTER(_mScore);
	setValues(_mScore, mScore, 1e100);

	for (int i = 1; i <= n1; i++)
		for (int j = 1; j <= n2; j++)
		{
			mScore[i-1 + n1 * (j-1)] = _afp_score(D1, D2, nD1, nD2, afpLength, i, j);
		}
	
	UNPROTECT(3);
	return (_mScore);
}

inline double _afp_dist(double *D1, double *D2, int nD1, int nD2, int afpLength, int i1, int i2, int j1, int j2)
{
	double d = 0;
	for (int k = 0; k < afpLength; k++)
	{
		d += R_pow_di(D1[i1+k-1 + nD1 * (i2+k-1)] - D2[j1+k-1 + nD2 * (j2+k-1)], 2);
	}
	d = sqrt(d/afpLength);
	return (d);
}

SEXP AFP_Distance(SEXP _D1, SEXP _D2, SEXP _nNodes, SEXP _nAFP, SEXP _afpLength)
{
	PROTECT(_D1 = AS_NUMERIC(_D1));
	PROTECT(_D2 = AS_NUMERIC(_D2));
	double *D1 = NUMERIC_POINTER(_D1);
	double *D2 = NUMERIC_POINTER(_D2);

	int nNodes = INTEGER_POINTER(AS_INTEGER(_nNodes))[0];
	int nAFP = INTEGER_POINTER(AS_INTEGER(_nAFP))[0];
	int afpLength = INTEGER_POINTER(AS_INTEGER(_afpLength))[0];

	SEXP _mDist;
	PROTECT(_mDist = NEW_NUMERIC(nAFP * nAFP * afpLength * (nNodes-afpLength)));
	setDim4(_mDist, nAFP, nAFP, afpLength, nNodes-afpLength);
	double *mDist = NUMERIC_POINTER(_mDist);
	setValues(_mDist, mDist, 1e100);

	int nD1 = nNodes, nD2 = nAFP+afpLength-1;

	int i1Max, iMax, j1Max;
	i1Max = nNodes - afpLength;
	for (int i1 = 1; i1 <= i1Max; i1++)
	{
		iMax = afpLength < i1? afpLength : i1;
		for (int i = 1; i <= iMax; i++)
		{
			j1Max = nAFP - i;
			for (int j1 = 1; j1 <= j1Max; j1++)
			{
				for (int j2 = j1+i; j2 <= nAFP; j2++)
				{
					mDist[j1-1 + nAFP * (j2-1 + nAFP * (i-1 + afpLength * (i1-1)))] = _afp_dist(D1, D2, nD1, nD2, afpLength, i1-i+1, i1+1, j1, j2);
				}
			}
		}
	}

	UNPROTECT(3);
	return (_mDist);
}
