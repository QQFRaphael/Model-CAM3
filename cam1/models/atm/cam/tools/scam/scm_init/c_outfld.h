#ifndef C_OUTFLD_H
#define C_OUTFLD_H

#include "realtype.h"
#include "fortran.h"
#include "ipc.h"

void    get_field( int id, real_t* outdata );
int     get_field_info( int id, field_info* fi );
void    reset_field( int id );
void    set_field( int id,  const real_t* indata );


#endif  /* C_OUTFLD_H */
