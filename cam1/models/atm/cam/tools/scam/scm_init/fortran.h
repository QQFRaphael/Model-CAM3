/*---------------------------------------------------------------------
 * DESCRIPTION:  Causes the preprocessor to modify all of the 
 *               FORTRAN function names in the C code to match
 *               the name assigned by the FORTRAN compiler, e.g.,
 *               append a trailing underscore
 */

#ifndef _Fortran_h
#define _Fortran_h

#if ( defined SUNOS ) || ( defined IRIX64 ) || ( defined OSF1 ) || ( defined LINUX )
#define FORTRAN_UNDERSCORE
#endif


#ifdef FORTRAN_UNDERSCORE
#define addfield   addfield_
#define init_c_outfld init_c_outfld_
#define bldfld bldfld_
#define get_levels get_levels_
#define init_model init_model_
#define c_outfld c_outfld_
#define stepon stepon_
/* following are for configuration tests */
#define run_tests run_tests_
#define test_int test_int_
#define test_double test_double_
#define test_int77 test_int77_
#define test_double77 test_double77_
#define test_int90 test_int90_
#define test_double90 test_double90_
#endif

#if ( defined LINUX )
#if ( defined PGF90 )
#define main MAIN_
#else
#define main MAIN__
#endif
#endif

#endif                          /* _Fortran_h */
/* DON'T ADD ANYTHING AFTER THIS #endif */








