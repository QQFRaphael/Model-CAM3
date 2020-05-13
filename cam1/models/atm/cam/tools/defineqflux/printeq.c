/*
** Jim Rosinski
** Student ID 006688918
** CS 5606
** Homework 6
**
** Function to print system of equations.
*/

#include <stdio.h>
#include "defineqflux.h"

void printeq (double *a[],            // A of matrix equation Ax=b
	      const int n,	     // rank of A
	      const double b[])	     // b vector of Ax=b
{
  int i;  // row index
  int k;  // column index

  for (i = 0; i < n; i++) {
    for (k = 0; k < n; k++)
      printf ("%9g ", a[i][k]);
    printf ("     %9g\n", b[i]);
  }
  printf ("\n");
}

void print_matrix (double **matrix,  // matrix to be printed
		   const int n)	    // rank of matrix
{
  int i;  // row index
  int k;  // column index

  for (i = 0; i < n; i++) {
    for (k = 0; k < n; k++)
      printf ("%9g ", matrix[i][k]);
    printf ("\n");
  }
  printf ("\n");
}

void print_vector (double *x,
		   const int n)
{
  int k;

  for (k = 0; k < n; k++)
    printf ("%lg ", x[k]);
  printf ("\n");
}
