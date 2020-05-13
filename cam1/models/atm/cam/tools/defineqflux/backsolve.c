/*
** Jim Rosinski
** Student ID 006688918
** CS 5606
** Homework 6
**
** Function to backsolve an upper triangular system (Ax=b)
** Follows discussion from class.
**
** Return value: 0 (success) or -1 (failure)
*/

#include <math.h>
#include "defineqflux.h"

int backsolve (double **a,        // matrix (input)
	       const int n,      // rank of matrix (input)      
	       const double b[],  // RHS of equation vector (input)
	       double x[])        // solution vector (output)
{
  int j;      // row loop index
  int k;      // loop index used for both rows and columns

  double sum;  // Sum of A(k,j)*x(j)

  // Loop over rows starting at the bottom of the triangular system

  for (k = n-1; k >= 0; k--) {
    if (a[k][k] == 0.) {
      printf ("Zero diagonal element %d\n", k);
      return -1;
    }

    // Note that for the 1st iteration (k=n-1) the sum is zero

    sum = 0.;
    for (j = k+1; j < n; j++)
      sum += a[k][j] * x[j];

    x[k] = (1. / a[k][k]) * (b[k] - sum);
  }

  return 0;
}
