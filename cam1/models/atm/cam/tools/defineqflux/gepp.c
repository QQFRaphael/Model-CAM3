/*
** Jim Rosinski
** Student ID 006688918
** CS 5606
** Homework 6
**
** Function to perform Gaussian elimination with partial pivoting.
** Follows discussion from class.  Note that indexing is zero-based.
**
** Return value: 0 (success) or -1 (failure)
*/

#include <math.h>
#include "defineqflux.h"

#define SWAP(a,b,tmp) (tmp = a); (a = b); (b = tmp);      // interchange elements

int gepp (double **a,     // matrix (input/output)
	  double **L,     // saved a(i,k)/a(k,k) values (input/output: input matrix ASSUMED zero)
	  const int n,    // rank of matrix (input)
	  double b[],     // B vector of Ax = b (input/output)
	  int p[])        // permutation vector P (output)
{
  int i;                  // row loop index
  int j;                  // column loop index
  int k;                  // loopi index used for both row and colun
  int r;                  // row loop index for pivoting
  int rmax;               // row index with max pivot element (the pivot row)
  int pswap;              // intermediate for swapping P vector elements

  double z;                // A(i,k)/A(k,k) (Fortran notation)
  double bswap;            // intermediate for swapping B vector elements

  double *tmp;         // row of A (struct contains pointer to row elements)

  // Initialize permuation vector

  for (k = 0; k < n; k++)
    p[k] = k;

  // Gaussian elimination loop

  for (k = 0; k < n-1; k++) {

    // Find pivot row

    rmax = k;
    for (r = k+1; r < n; r++) {
      if (fabs (a[r][k]) > fabs (a[rmax][k]))
	rmax = r;
    }

    if (a[rmax][k] == 0.) {
      printf ("Row %d col %d of A became zero => singular\n", p[r], k);
      return -1;
    }

    // Do the pivoting (if required) of matrix A rows and vector B. Update pivot 
    // vector P. The swapping of rows of A only involves a pointer switch.

    if (rmax != k) {
      printf ("Permuting rows %d and %d\n", k, rmax);
      SWAP (a[k], a[rmax], tmp);
      SWAP (L[k], L[rmax], tmp);
      SWAP (b[k], b[rmax], bswap);
      SWAP (p[k], p[rmax], pswap);
    }

#ifdef DEBUG
    printf ("System of equations after pivoting k=%d:\n", k);
    printeq (a, n, b);
    printf ("L matrix:\n");
    print_matrix (L, n);
#endif

    // Redefine row i of the matrix.  Explicitly zeroing matrix elements
    // which would come out that way anyway is helpful if printing
    // the intermediate system of equations (DEBUG defined).

    for (i = k+1; i < n; i++) {
      z = a[i][k] / a[k][k];
      L[i][k] = z;
#ifdef DEBUG
      for (j = 0; j < k+1; j++)
	a[i][j] = 0.;
#endif
      for (j = k+1; j < n; j++)
	a[i][j] -= z*a[k][j];

      b[i] -= z*b[k];
    }
#ifdef DEBUG
    printf ("System of equations after G.E. k=%d:\n", k);
    printeq (a, n, b);
    printf ("L matrix:\n");
    print_matrix (L, n);
#endif
  }

  if (a[n-1][n-1] == 0.) {
    printf ("Diagonal element %d of A is zero => singular\n", p[n-1]);
    return -1;
  }
  return 0;
}
