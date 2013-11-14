
#include "common.h"

/* conditional mutual information, to be used in C code. */
double c_cmi(int *n, int lx, int ly, int lz) {

int xi, yi, zi, sum; 
double res = 0;

  for (zi = 0; zi < lz; zi++)
    for (yi = 0; yi < ly; yi++)
      for (xi = 0; xi < lx; xi++)
        res += MI_PART(
          n[xi + yi*(lx+1) + zi*(ly+1)*(lx+1)], /*xyz joint*/
          n[xi + ly*(lx+1) + zi*(ly+1)*(lx+1)], /*xz margin*/
          n[lx + yi*(lx+1) + zi*(ly+1)*(lx+1)], /*yz margin*/
          n[lx + ly*(lx+1) + zi*(ly+1)*(lx+1)]); /*z margin*/

  /*trick to handle both 2d and 3d tables*/
  if (lz == 1)
    sum = n[lx + ly*(lx+1)];
  else
    sum = n[lx + ly*(lx+1) + lz*(ly+1)*(lx+1)];

  return res / sum;

}/*C_CMI*/

/* unconditional mutual information, to be used for the asymptotic test. */
SEXP mi(SEXP x, SEXP y, SEXP gsquare, SEXP df_adjust) {

int lx = NLEVELS(x), ly = NLEVELS(y);
int num = LENGTH(x), adj = LOGICAL(df_adjust)[0];
int *xx = INTEGER(x), *yy = INTEGER(y), *n;
double *res = NULL;
SEXP result;

  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);

  /* build the contingency table. */
  n = table_2d(xx, lx, yy, ly, num);

  /* compute statistic and df. */
  res[0] = c_cmi(n, lx, ly, 1);
  res[1] = c_df(n, lx, ly, 1, adj);

  /* rescale to match the G^2 test. */
  if (isTRUE(gsquare))
    res[0] *= 2 * num;

  UNPROTECT(1);

  return result;

}/*MI*/

/* conditional mutual information, to be used for the asymptotic test. */
SEXP cmi(SEXP x, SEXP y, SEXP z, SEXP gsquare, SEXP df_adjust) {

int lx = NLEVELS(x), ly = NLEVELS(y), lz = NLEVELS(z);
int num = LENGTH(x), adj = LOGICAL(df_adjust)[0];
int *xx = INTEGER(x), *yy = INTEGER(y), *zz = INTEGER(z), *n;
double *res = NULL;
SEXP result;

  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);

  /* build the contingency table. */
  n = table_3d(xx, lx, yy, ly, zz, lz, num);

  /* compute statistic and df. */
  res[0] = c_cmi(n, lx, ly, lz);
  res[1] = c_df(n, lx, ly, lz, adj);

  /* rescale to match the G^2 test. */
  if (isTRUE(gsquare))
    res[0] *= 2 * num;

  UNPROTECT(1);

  return result;

}/*CMI*/

/* unconditional Gaussian mutual information, to be used in C code. */
double c_mig(double *xx, double *yy, int *num) {

double cor = c_fast_cor(xx, yy, num);

  return - 0.5 * log(1 - cor * cor);

}/*C_MIG*/

/* unconditional Gaussian mutual information, to be used in the asymptotic test. */
SEXP mig(SEXP x, SEXP y, SEXP length) {

double *xx = REAL(x), *yy = REAL(y);
int *num = INTEGER(length);
SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));
  NUM(result) = c_mig(xx, yy, num);
  UNPROTECT(1);

  return result;

}/*MIG*/
