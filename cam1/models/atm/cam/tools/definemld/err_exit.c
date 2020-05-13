/*
** Generic "error exit" program to print a message and exit with a non-zero
** status.
*/
 
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
 
void err_exit (const char *fmt, ...)
{
  va_list args;
 
  va_start (args, fmt);
  if (fmt != NULL) {
    vfprintf (stderr, fmt, args);
  }
  va_end (args);
  abort ();
  exit (1);
}
