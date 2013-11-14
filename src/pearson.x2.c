
#include "common.h"

/* conditional Pearson's x2, to be used in C code. */
double c_cx2(int *n, int lx, int ly, int lz) {

int xi, yi, zi, sum; 
double res = 0;

  for (zi = 0; zi < lz; zi++)
    for (yi = 0; yi < ly; yi++)
      for (xi = 0; xi < lx; xi++)
        res += X2_PART(
          n[xi + yi*(lx+1) + zi*(ly+1)*(lx+1)], /*xyz joint*/
          n[xi + ly*(lx+1) + zi*(ly+1)*(lx+1)], /*xz margin*/
          n[lx + yi*(lx+1) + zi*(ly+1)*(lx+1)], /*yz margin*/
          n[lx + ly*(lx+1) + zi*(ly+1)*(lx+1)]); /*z margin*/

  return res;

}/*C_CX2*/

SEXP x2 (SEXP x, SEXP y, SEXP df_adjust) {

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
  res[0] = c_cx2(n, lx, ly, 1);
  res[1] = c_df(n, lx, ly, 1, adj);

  UNPROTECT(1);

  return result;

}/*X2*/

SEXP cx2 (SEXP x, SEXP y, SEXP z, SEXP df_adjust) {

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
  res[0] = c_cx2(n, lx, ly, lz);
  res[1] = c_df(n, lx, ly, lz, adj);

  UNPROTECT(1);

  return result;

}/*CX2*/

