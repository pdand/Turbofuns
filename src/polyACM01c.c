#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

void F77_NAME(polyACM_f) (int *ncase, int *nvar, int *IAdjust, int *NCore, int *iRaw, double*xThreshold, double *xPoly, int *IError, double *xACM);

void F77_NAME(polyR_f) (int *ncase, int *nvar, int *IAdjust, int *NCore, int *iRaw, double*xThreshold, double *xPoly, int *IError);


  extern SEXP c_polyACM_f(SEXP ncase, SEXP nvar, SEXP IAdjust, SEXP NCore, SEXP iRaw){

   SEXP oxThreshold = PROTECT(allocMatrix(REALSXP, 11, asInteger(nvar)));
   SEXP oxPoly = PROTECT(allocMatrix(REALSXP, asInteger(nvar),asInteger(nvar)));
   SEXP oIError = PROTECT(allocMatrix(INTSXP, 2, asInteger(nvar)*(asInteger(nvar)-1)/2));
   SEXP oxACM = PROTECT(allocMatrix(REALSXP, asInteger(nvar)*(asInteger(nvar)-1)/2,asInteger(nvar)*(asInteger(nvar)-1)/2));
   F77_CALL(polyACM_f) (INTEGER(ncase), INTEGER(nvar), INTEGER(IAdjust), INTEGER(NCore), INTEGER(iRaw), REAL(oxThreshold), REAL(oxPoly), INTEGER(oIError),REAL(oxACM));
   SEXP OS = PROTECT(allocVector(VECSXP, 4));
   SET_VECTOR_ELT(OS,0,oxThreshold);
   SET_VECTOR_ELT(OS,1,oxPoly);
   SET_VECTOR_ELT(OS,2,oIError);
   SET_VECTOR_ELT(OS,3,oxACM);
   UNPROTECT(5);
   return OS;
 }


   extern SEXP c_polyR_f(SEXP ncase, SEXP nvar, SEXP IAdjust, SEXP NCore, SEXP iRaw){

   SEXP oxThreshold = PROTECT(allocMatrix(REALSXP, 11, asInteger(nvar)));
   SEXP oxPoly = PROTECT(allocMatrix(REALSXP, asInteger(nvar),asInteger(nvar)));
   SEXP oIError = PROTECT(allocMatrix(INTSXP, 2, asInteger(nvar)*(asInteger(nvar)-1)/2));
   F77_CALL(polyR_f) (INTEGER(ncase), INTEGER(nvar), INTEGER(IAdjust), INTEGER(NCore), INTEGER(iRaw), REAL(oxThreshold), REAL(oxPoly), INTEGER(oIError));
   SEXP OS = PROTECT(allocVector(VECSXP, 3));
   SET_VECTOR_ELT(OS,0,oxThreshold);
   SET_VECTOR_ELT(OS,1,oxPoly);
   SET_VECTOR_ELT(OS,2,oIError);
   UNPROTECT(4);
   return OS;
 }


 static const R_CallMethodDef CallEntries[] = {
   {"c_polyACM_f",   (DL_FUNC) &c_polyACM_f,   5},
   {"c_polyR_f",   (DL_FUNC) &c_polyR_f,   5},
   {NULL,         NULL,                0}
 };



void R_init_Turbofuns(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  }

/*
void R_init_TestMFun(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
*/
