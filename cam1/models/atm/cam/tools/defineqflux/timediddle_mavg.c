#include "defineqflux.h"

void timediddle_mavg (double qflux[12][PLAT][PLON])
{
  extern const int daysmo[12];
    
  int i, j;        // lon, lat indices
  int m, mm, mp;   // month indices 
  int n;           // month index
  int p[12];       // pivot vector

  double ddaysmo[12];
  double **A;      // 2-d matrix: define as pointers to accomodate gepp
  double **L;      // 2-d lower triangular matrix for gepp
  double b[12];    // Ax = b matrix equation.  b is known (monthly mean qflux)
  double x[12];    // Ax = b matrix equation.  x is the solution vector (time diddled)

  if ((A = (double **) malloc (12 * sizeof (double *))) == 0)
    err_exit ("malloc failure for A\n");

  if ((L = (double **) malloc (12 * sizeof (double *))) == 0)
    err_exit ("malloc failure for L\n");

  for (m = 0; m < 12; m++) {

    if ((A[m] = (double *) malloc (12 * sizeof (double))) == 0)
      err_exit ("malloc failure for row %d of A\n", m);

    if ((L[m] = (double *) malloc (12 * sizeof (double))) == 0)
      err_exit ("malloc failure for row %d of A\n", m);

    ddaysmo[m] = daysmo[m];
  }

  for (j = 0; j < PLAT; j++) {
    for (i = 0; i < PLON; i++) {
      for (m = 0; m < 12; m++) {

	for (n = 0; n < 12; n++) {
	  A[m][n] = 0.;
	  L[m][n] = 0.;
	}
	mm = m == 0 ? 11 : m - 1;
	mp = (m + 1) % 12;
	A[m][m]   = 1. - 0.25*(ddaysmo[m]/(ddaysmo[m] + ddaysmo[mm]) + 
			       ddaysmo[m]/(ddaysmo[m] + ddaysmo[mp]));
	A[m][mm] = 0.25*(ddaysmo[m]/(ddaysmo[m]+ddaysmo[mm]));
	A[m][mp] = 0.25*(ddaysmo[m]/(ddaysmo[m]+ddaysmo[mp]));
	b[m] = qflux[m][j][i];
      }

      // First do the Gaussian elimination with partial pivoting
  
      if (gepp (A, L, 12, b, p) < 0)
	err_exit ("gepp failure\n");

      // Now backsolve to get the solution vector, X

      if (backsolve (A, 12, b, x) < 0)
	err_exit ("backsolve failure\n");

      // Store the solution vector back in qflux

      for (m = 0; m < 12; m++) {
	qflux[m][j][i] = x[m];
      }
    }
  }
  
  return;
}
