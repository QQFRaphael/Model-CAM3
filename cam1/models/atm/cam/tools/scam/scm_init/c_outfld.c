#include <params.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "ipc.h"
#include "fortran.h"
#include "max.h"

#if ( defined FORTRANCAPS )
#define c_outfld C_OUTFLD
#elif ( defined FORTRANUNDERSCORE )
#define c_outfld c_outfld_
#elif ( defined FORTRANDOUBLEUNDERSCORE )
#define c_outfld c_outfld__
#endif

#define MAX_HASH (94 * 4)       /* highest value returned by 4 ascii chars */
#define MAX_DEPTH (100)          /* max number of fields that can share the same hash value */


#if (defined USE_4BYTE_REAL )
typedef float   ccm_real_t;
#else
typedef double   ccm_real_t;
#endif

/*
 * field struct to hold all information related to storage of field data
 */
static struct field {
    int    is_averaged;
      /* is the field averaged */
    int    is_modifiable;
      /* is the field modifiable */
    int    is_shown;
      /* is the field shown in the user interface */
    char   name[MAX_FIELDNAME_LEN];
      /* name of the field */
    char   longname[MAX_FIELDNAME_LEN];
      /* full name of the field */
    ccm_real_t min;
      /* expected minimum range of field data */
    ccm_real_t max;
      /* expected maximum range of field data */    
    ccm_real_t mult;
      /* multiplication factor for displaying data in useful units */
    char   units[20];
      /* units to use when displaying data after multiplier is applied */
    char   std_units[20];
      /* standard si units used in the model */    
    int    size;
      /* size of field in bytes */
    int    numouts;
      /* number of times c_outfld has been called since last get_field call */
    ccm_real_t data[MAX_LEVELS];
      /* buffer to hold data from c_outfld call */
    ccm_real_t* ptr;
      /* pointers to modifiable fields */
};

typedef struct field field;

static int field_count = 0;     /* number of fields created by bldfld */



/*
 * hash table for fast retrieval of fields by name
 */
static struct field_hash_table {
    field* fields[MAX_HASH][MAX_DEPTH];
    int    depth[MAX_HASH];
} ftable;

/*
 *  array for fast retrieval of fields by id
 */ 
static struct field_array {
    field* fields[MAX_FIELDS];
} farray;



/*
 *   hash function for computing hash value of field names
 */ 
int
hash_func( const char* name )
{
    int hash = 0;
    int i = 0;

      /* hash on first 4 characters of the name */
    for ( ; i<4; i++ ){
        if ( name[i] == 0 )
            break;
        hash += (int)name[i] - 32; /* ascii values of letters begins at 32 */ 
    }
    return hash;
}


/*
 *  retrieves a field by hash lookup of field name
 */
field*
lookup_field( const char* name )
{
    int i;
    
    int hash = hash_func( name );
    for ( i=0; i<ftable.depth[hash]; i++) {
        field* f = ftable.fields[hash][i];
        if ( ! strcmp( name, f->name ) )
            return f;
    }
      /* didn't find it */
    return (field*) 0;
}


/*
 *  retrieves a field by array indexing of id
 */
field*
lookup_field_by_id( int id )
{
    return farray.fields[id];
}

    
/*
 *   insert a field into both the hash table (by name) and the array (by id)
 */
void
insert_field( const char* name, field* f )
{
    int hash = hash_func( name );
    if ( ftable.depth[hash] == MAX_DEPTH ) {
        fprintf( stderr, "ERROR: HASH TABLE FULL - you must increase" 
                 " the value of MAX_DEPTH in c_outfld.c" );
        exit( 1 );
    }
    ftable.fields[hash][ftable.depth[hash]++] = f;
      /*
       * this should be modified to dynamically increase the size of the hash table
       * or the hash table should use a linked list, etc.
       */
    farray.fields[field_count] = f;
    field_count++;
}


/*
 *  add a field to the fields hashtable and array, ignoring duplicates
 */
void
addfield( char* f_name, char* longname, char* units,
          char* std_units, ccm_real_t* mult, ccm_real_t* min,         
          ccm_real_t*  max, int* is_shown, int* is_modifiable,    
          int* is_averaged, int*  size, ccm_real_t* dummy,
          int strlen1, int strlen2,
          int strlen3, int strlen4 )
{
    
    int i;
    static int field_idx = 0;
    char c_name[MAX_FIELDNAME_LEN];
    field* f;
    
      /*
       * determine actual length of string: strlen includes blanks
       */
    while ( isspace(f_name[strlen1-1]) )
        strlen1--;

    strncpy( c_name, f_name, strlen1 );
    c_name[strlen1] = 0;
      /*
       * warn if addfield called twice with same field
       */
    if ( lookup_field( c_name ) != 0 ) {
        fprintf( stderr, "WARNING: addfield called twice for \"%s\"\n", c_name ); 
    }

    /*    printf("name of field %s and size of field is %i\n",c_name,*size);*/
    f = (field*) malloc( sizeof( field ) );
    /*    printf("address of f is %x size of field is %i\n",f, sizeof( field )); */
    if ( ! f ) {
        fprintf( stderr,"ERROR: addfield: out of memory" );
        exit( -1 );
    }
                 
    strcpy( f->name, c_name );
    strncpy( f->longname, longname, strlen2 );
    f->longname[strlen2] = 0;
    strncpy( f->units, units, strlen3 );
    f->units[strlen3] = 0;
    strncpy( f->std_units, std_units, strlen4 );
    f->std_units[strlen4] = 0;
    f->mult = *mult;
    f->min = *min;
    f->max = *max;
    f->is_shown = *is_shown;
    f->is_averaged = *is_averaged;
    f->is_modifiable = *is_modifiable;
    f->size = *size;
    for ( i=0; i<*size; i++ )
        f->data[i] = 0;
    f->ptr = (ccm_real_t*)0;
    f->numouts = 0;
    insert_field( f->name, f ); /* add the field to the hash table */
}


/***************************************************************************
 *  get_field:  finds field by id, copies data from the field data buffer 
 * 
 *  inputs:
 *      id  -    id of desired field
 *  outputs:
 *      data  -    data copied from the field data buffer
 ***************************************************************************/

void
get_field( int id,  real_t* data )
{
    int i;
    
    field* f = lookup_field_by_id( id ); 

    if ( f == 0 ) {
        fprintf(stderr,"Error: get_field: field id number %d not found\n", id );
        exit ( 1 );
    }
      /* check that modifiable fields pointer is initialized - this should be done by
       * making an c_outfld call in intht.F during the model initialization
       * process
       */
    if ( f->is_modifiable )
        if ( !f->ptr ) {
            fprintf(stderr,"Error: get_field() - modifiable field %s data pointer is not initialized!\nMake sure that intht.F contains an c_outfld call for this field.\n", f->name );
            exit ( 1 );
        }

      /*
       * if the field is averaged, divide the value by the number of times
       * c_outfld has been called since the last time the data was retrieved
       */
    
    if ( f->is_averaged ) {
        if ( f->numouts == 0 )
            for ( i=0; i<f->size; i++ ) 
                data[i] = f->data[i];
        else
            for ( i=0; i<f->size; i++ ) 
                data[i] = f->data[i]/f->numouts;  
    }
    else {
        if ( sizeof(*data) == sizeof(ccm_real_t) )
            memcpy( data, f->data, f->size*sizeof(ccm_real_t) );
        else
            for ( i=0; i<f->size; i++ ) 
                data[i] = f->data[i];
    }
}

 

/***************************************************************************
 *  get_field_info:  retrieves field information 
 * 
 *  inputs:
 *     id - the id of the target field
 *  returns: 0 if field id exists, -1 otherwise
 ***************************************************************************/

int
get_field_info( int id, field_info* fi )
{
    field* f;

    if ( id >= field_count ) {
        fi->size = 0;           
        return -1;
    }
    
    f = lookup_field_by_id( id );
    
    strcpy ( fi->name, f->name );
    
    strcpy ( fi->longname, f->longname );
    strcpy ( fi->units, f->units );
    strcpy ( fi->std_units, f->std_units );
    fi->mult = f->mult;
    fi->min = f->min;
    fi->max = f->max;
    fi->is_shown = f->is_shown;
    fi->is_modifiable = f->is_modifiable;
    fi->is_averaged = f->is_averaged;
    fi->size = f->size;
    return 0;
}


/***************************************************************************
 *  reset_field:  find field by name, reset the numouts, zero out
 *                 the data buffer.
 * 
 *  inputs:
 *      id - the id of the field
 *
 ***************************************************************************/

void
reset_field( int id )
{
    field* f = lookup_field_by_id( id );
    f->numouts = 0;
    memset( f->data, 0, f->size*sizeof(ccm_real_t) );
}


/***************************************************************************
 *  set_field:  find field by id, copie data to the field buffer, pointer 
 * 
 *  inputs:
 *      id  - id of the field
 *      data  - data to be copied 
 ***************************************************************************/

void
set_field( int id,  const real_t* data )
{
    int i;
    field* f = lookup_field_by_id( id );

    if ( f == 0 ) {
        fprintf(stderr,"Error: get_field: field number %d not found\n", id );
        exit ( 1 );
    }

    /*
     * need to copy the data to both the field data buffer
     * and the memory location in the model (if modifiable field)
     */
    
      /* check that pointer is initialized */
    if ( f->is_modifiable )
        if ( !f->ptr ) {
            fprintf(stderr,"Error: set_field() - modifiable field %s data pointer is not initialized!\nMake sure that intht.F contains an c_outfld call for this field.\n", f->name );
            exit ( 1 );
        }
        
    if ( sizeof(*data) == sizeof(ccm_real_t) ) {
        memcpy( f->data, data, f->size*sizeof(ccm_real_t) );
        if ( f->is_modifiable )
	  {
            memcpy( f->ptr, data, f->size*sizeof(ccm_real_t) );
	    //	    	    printf("%f ************\n",data[0]); 
	    //	    	    printf("%d ****fptr********\n",f->ptr);
	  }
    } 
    else {
        for ( i=0; i<f->size; i++ ) {
            f->data[i] = data[i];
            if ( f->is_modifiable )
                f->ptr[i] = data[i];
        }
    }
    f->numouts = 1;
}


/*****************************************************************************
 *
 *  c_outfld:  finds fields by name, copies data into field data buffer
 *
 *  INPUTS:
 *      name :   The name of the field to be output;
                    called by fortran, so no null terminator
 *      data :     data to be copied to the buffer
 *      plon, lat:  Not used, kept for compatibility with
 *                         original CCM procedure call.
 *      namelen:   The length of the name.In most fortran-C interfaces,
 *                 the lengths of char strings passed as parameters are
 *                 invisibly added to the end of the arg list.
 *
 ******************************************************************************/ 
void
c_outfld( char* f_name, ccm_real_t* data, int* plon, int* lat, int namelen )
{
    int    i,k;
    field* f;
    char c_name[MAX_FIELDNAME_LEN];
    
      /* ignore blanks */
    while ( isspace(f_name[namelen-1]) )
        namelen--;
    strncpy( c_name, f_name, namelen );
    c_name[namelen] = 0;

      /* perform hash table lookup by name */
    f = lookup_field( c_name );
    if ( f == 0 ) {
        fprintf(stderr,"Warning: c_outfld: field \"%s\" not found\n", c_name );
        return;
    }
      /*
       * if this is a modifiable field then save the pointer (needed for setting data)
       */
    if ( f->is_modifiable && !f->ptr )
        f->ptr = data;
      /*
       * Copy the data to the buffer 
       */
    f->numouts++;
    if ( f->is_averaged ) {
        for ( k=0; k<f->size; k++ )
            f->data[k] += data[k];
    }
    else {
        if ( sizeof(*data) == sizeof(ccm_real_t) )
            memcpy( f->data, data, f->size*sizeof(ccm_real_t) );
        else
            for ( k=0; k<f->size; k++ ) 
                f->data[k] = data[k];
    }
    return;
}

/*
 * init_c_outfld: initialize c_outfld data structures
 */
void
init_c_outfld()
{
    int i;
    
    for ( i=0; i<MAX_HASH; i++ )
        ftable.depth[i] = 0;
      /* free up memory used by previous addfield calls */
    for ( i=0; i<field_count; i++ )
        free ( (char*)lookup_field_by_id( i ));
    field_count = 0;
}




