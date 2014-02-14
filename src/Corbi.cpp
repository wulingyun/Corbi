#include "Corbi.h"

const R_CallMethodDef callMethods[] = {
  {"RMultiHyper", (DL_FUNC) &RMultiHyper, 3},
  {"PA_Scores", (DL_FUNC) &PA_Scores, 5},
  {"NQ_ShortestDistances", (DL_FUNC) &NQ_ShortestDistances, 2},
  {"NE_Depths", (DL_FUNC) &NE_Depths, 2},
  {NULL, NULL, 0}
};

void R_init_Corbi(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
