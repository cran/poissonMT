#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(outerr)(void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"outerr", (DL_FUNC) &F77_NAME(outerr), 4},
    {NULL, NULL, 0}
};

void R_init_poissonMT(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
