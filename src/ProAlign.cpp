#include "Corbi.h"
#include <Rmath.h>

double calc_dist(double *D1, double *D2, int i1, int i2, int j1, int j2, int n1, int n2, int afpLength)
{
	double d = 0;
	for (int i = 0; i < afpLength; i++)
	{
		d += R_pow_di(D1[i1+i-1 + n1 * (i2+i-1)] - D2[j1+i-1 + n2 * (j2+i-1)], 2);
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
					mDist[j1-1 + nAFP * (j2-1 + nAFP * (i-1 + afpLength * (i1-1)))] = calc_dist(D1, D2, i1-i+1, i1+1, j1, j2, nNodes, nAFP+afpLength-1, afpLength);
				}
			}
		}
	}

	UNPROTECT(3);
	return (_mDist);
}
