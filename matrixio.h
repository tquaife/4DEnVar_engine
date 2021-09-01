#ifndef MATRIXIO_H
#define MATRIXIO_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

#define MATRIXIO_MAX_LINE_LEN 50000


void	alloc_error( char *caller_name );

/*
** Long int functions
*/

char	get_ldElement( char *Line, long *Element_p );
long 	*read_ascii_ldMatrix(  long *X_size, long *Y_size );
long 	*fread_ascii_ldMatrix(  long *X_size, long *Y_size, FILE * );
void	print_ldArray( long *Array_p, long size );
void	print_ldMatrix( long *Matrix_p, long x_size, long y_size );

/*
** Double functions
*/

char	get_ffElement( char *Line, double *Element_p );
double 	*read_ascii_ffMatrix(  long *X_size, long *Y_size );
double 	*fread_ascii_ffMatrix(  long *X_size, long *Y_size, FILE *stream );
void	print_ffArray( double *Array_p, long size );
void	print_ffMatrix( double *Matrix_p, long x_size, long y_size );

/*
** float functions
*/

char	get_fElement( char *Line, float *Element_p );
float 	*read_ascii_fMatrix(  long *X_size, long *Y_size );
float 	*fread_ascii_fMatrix(  long *X_size, long *Y_size, FILE *stream );
void	print_fArray( float *Array_p, long size );
void	print_fMatrix( float *Matrix_p, long x_size, long y_size );
void	flip_fMatrix( float *Matrix_P, long x_size, long y_size );


#endif


