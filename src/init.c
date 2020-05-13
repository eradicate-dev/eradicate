#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP getDetVecs(SEXP y_arr, SEXP mp_arr, SEXP J_i, SEXP tin, SEXP K_) ;
extern SEXP getSingleDetVec(SEXP y_, SEXP mp_, SEXP K_);
static const R_CallMethodDef CallEntries[] = {
    {"getDetVecs",      (DL_FUNC) &getDetVecs,       5},
    {"getSingleDetVec", (DL_FUNC) &getSingleDetVec,  3},
    {NULL, NULL, 0}
};

void R_init_eradicate(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
