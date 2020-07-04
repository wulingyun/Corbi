#include "Corbi.h"

const R_CallMethodDef callMethods[] = {
  {"RMultiHyper", (DL_FUNC) &RMultiHyper, 3},
  {"PMultiHyper", (DL_FUNC) &PMultiHyper, 4},
  {"PMultiNom", (DL_FUNC) &PMultiNom, 4},
  {"PA_Scores", (DL_FUNC) &PA_Scores, 5},
  {"NQ_ShortestDistances", (DL_FUNC) &NQ_ShortestDistances, 3},
  {"BS_GetSubnets", (DL_FUNC) &BS_GetSubnets, 3},
  {"BS_ExtendSubnets", (DL_FUNC) &BS_ExtendSubnets, 3},
  {"ND_RatioDistribution", (DL_FUNC) &ND_RatioDistribution, 2},
  {"ND_RatioDistributionParI", (DL_FUNC) &ND_RatioDistributionParI, 3},
  {"ND_ParMerge", (DL_FUNC) &ND_ParMerge, 4},
  {"ND_RatioDistributionAB", (DL_FUNC) &ND_RatioDistributionAB, 3},
  {"ND_RatioDistributionParAiB", (DL_FUNC) &ND_RatioDistributionParAiB, 3},
  {"ND_RatioDistribution2", (DL_FUNC) &ND_RatioDistribution2, 3},
  {"ND_DiffRatioNet", (DL_FUNC) &ND_DiffRatioNet, 2},
  {"ND_RatioVariance", (DL_FUNC) &ND_RatioVariance, 1},
  {"ND_RatioVarianceParI", (DL_FUNC) &ND_RatioVarianceParI, 1},
  {NULL, NULL, 0}
};

void R_init_Corbi(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
