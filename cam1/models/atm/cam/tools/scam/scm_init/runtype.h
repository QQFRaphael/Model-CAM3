
#ifndef RUNTYPE_H
#define RUNTYPE_H  /* this file is included in both FORTRAN and C code files */


#define ANAL    0
#define MODEL   1
#define IOP     2
#define LSMINI  3
#define LSMSRF  4
#define OZON    5
#define PRES    6
#define SST     7
#define ABSEMS  8
#define AEROPTICS 9
#define AERMASS  10
#define LSMPFT  11
#define USER    12
#define SIC     13               /* SIC **has** to be last */
#define NUM_INIT_FILES 14

#define MAX_RUNTYPE SIC + USER
#endif /*  RUNTYPE_H */







