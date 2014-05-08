#include "Corbi.h"

const R_CallMethodDef callMethods[] = {
  {"RMultiHyper", (DL_FUNC) &RMultiHyper, 3},
  {"PMultiHyper", (DL_FUNC) &PMultiHyper, 4},
  {"PA_Scores", (DL_FUNC) &PA_Scores, 5},
  {"NQ_ShortestDistances", (DL_FUNC) &NQ_ShortestDistances, 3},
  {"NE_GetDepths", (DL_FUNC) &NE_GetDepths, 4},
  {"NE_CountDepths", (DL_FUNC) &NE_CountDepths, 2},
  {NULL, NULL, 0}
};

void R_init_Corbi(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
