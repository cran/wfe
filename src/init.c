#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void CalDID(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void comp_OmegaHAC(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void comp_OmegaHC(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Demean(void *, void *, void *, void *, void *);
extern void DemeanDID(void *, void *, void *, void *, void *, void *, void *, void *);
extern void GenTime(void *, void *, void *, void *);
extern void GenWeightsDID(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GenWeightsFD(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GenWeightsMDID(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GenWeightsTime(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GenWeightsUnit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Index(void *, void *, void *, void *, void *);
extern void LamdaDID1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void LamdaDID2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void MDummy(void *, void *, void *, void *);
extern void OmegaDiDHAC(void *, void *, void *, void *, void *, void *, void *, void *);
extern void OmegaDiDHAC2(void *, void *, void *, void *, void *, void *, void *, void *);
extern void OmegaHatHAC(void *, void *, void *, void *, void *, void *, void *);
extern void OmegaHatHC(void *, void *, void *, void *, void *, void *, void *);
extern void ProjectionM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Transform(void *, void *, void *, void *, void *);
extern void TwayDemean(void *, void *, void *, void *, void *, void *, void *);
extern void VectorizeC(void *, void *, void *, void *, void *, void *, void *);
extern void WDemean(void *, void *, void *, void *, void *, void *);
extern void WWDemean(void *, void *, void *, void *, void *, void *);
extern void XWXiSum(void *, void *, void *, void *, void *, void *, void *);
extern void XXiSum(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"CalDID",         (DL_FUNC) &CalDID,         12},
    {"comp_OmegaHAC",  (DL_FUNC) &comp_OmegaHAC,   9},
    {"comp_OmegaHC",   (DL_FUNC) &comp_OmegaHC,    9},
    {"Demean",         (DL_FUNC) &Demean,          5},
    {"DemeanDID",      (DL_FUNC) &DemeanDID,       8},
    {"GenTime",        (DL_FUNC) &GenTime,         4},
    {"GenWeightsDID",  (DL_FUNC) &GenWeightsDID,  11},
    {"GenWeightsFD",   (DL_FUNC) &GenWeightsFD,   11},
    {"GenWeightsMDID", (DL_FUNC) &GenWeightsMDID, 13},
    {"GenWeightsTime", (DL_FUNC) &GenWeightsTime, 11},
    {"GenWeightsUnit", (DL_FUNC) &GenWeightsUnit, 11},
    {"Index",          (DL_FUNC) &Index,           5},
    {"LamdaDID1",      (DL_FUNC) &LamdaDID1,      14},
    {"LamdaDID2",      (DL_FUNC) &LamdaDID2,      14},
    {"MDummy",         (DL_FUNC) &MDummy,          4},
    {"OmegaDiDHAC",    (DL_FUNC) &OmegaDiDHAC,     8},
    {"OmegaDiDHAC2",   (DL_FUNC) &OmegaDiDHAC2,    8},
    {"OmegaHatHAC",    (DL_FUNC) &OmegaHatHAC,     7},
    {"OmegaHatHC",     (DL_FUNC) &OmegaHatHC,      7},
    {"ProjectionM",    (DL_FUNC) &ProjectionM,    17},
    {"Transform",      (DL_FUNC) &Transform,       5},
    {"TwayDemean",     (DL_FUNC) &TwayDemean,      7},
    {"VectorizeC",     (DL_FUNC) &VectorizeC,      7},
    {"WDemean",        (DL_FUNC) &WDemean,         6},
    {"WWDemean",       (DL_FUNC) &WWDemean,        6},
    {"XWXiSum",        (DL_FUNC) &XWXiSum,         7},
    {"XXiSum",         (DL_FUNC) &XXiSum,          6},
    {NULL, NULL, 0}
};

void R_init_wfe(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
