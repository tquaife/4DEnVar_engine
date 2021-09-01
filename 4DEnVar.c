/*

Requires my matrix i/o code and the GSL.

Tristan Quaife. 2021.
tquaife@gmail.com
*/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl_matrix.h>
#include<matrixio.h>
#include<4DEnVar_engine.h>

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/*
============================
Some basic utility functions
============================
*/


void check_dims( gsl_matrix * gmat1, gsl_matrix * gmat2 )
/*Check the dimensions of two GSL matrices.
Exit with failure if they are not equal.*/
{
	if( (gmat1->size1 != gmat2->size1) || (gmat1->size2 != gmat2->size2) ){
		fprintf( stderr, "%s: matrix dimensions do not match\n", __FILE__ );
		exit( EXIT_FAILURE );
	}
	return ;
}


gsl_matrix *load_ffMatrix_to_gsl_matrix( char *filename )
/*Wrapper function for fread_ascii_ffMatrix that handles
the file opening and exits with failure if open fails.
Returns the data as a GSL matrix.
*/
{
	FILE * fp ;
	double *ptr, tmp ;
	int x, y ;
	long xsize ; 
	long ysize ;
	
	gsl_matrix *gmat ;
	
	if( ( fp = fopen( filename, "r" ) ) == NULL ){
		fprintf( stderr, "%s: unable to open file %s\n", __FILE__, filename );
		exit( EXIT_FAILURE );
	}
	
	ptr=fread_ascii_ffMatrix( &xsize, &ysize, fp ) ;
    gmat=(gsl_matrix *)gsl_matrix_alloc ( (size_t)ysize, (size_t)xsize);

    for( y=0; y<ysize; y++ ){
        for( x=0; x<xsize; x++ ){
            tmp=*(ptr + y*xsize + x );
            gsl_matrix_set (gmat, y, x, tmp);
        }
    }

	fclose(fp);	
	free(ptr);
	return(gmat);
}

gsl_vector *load_ffMatrix_to_gsl_vector( char *filename )
/*Wrapper function for fread_ascii_ffMatrix that handles
the file opening and exits with failure if open fails.
Returns the data as a GSL VECTOR.
*/
{
	FILE * fp ;
	double *ptr, tmp ;
	int y ;
	long xsize ; 
	long ysize ;
	
	gsl_vector *gvec ;
	
	if( ( fp = fopen( filename, "r" ) ) == NULL ){
		fprintf( stderr, "%s: unable to open file %s\n", __FILE__, filename );
		exit( EXIT_FAILURE );
	}
	
	ptr=fread_ascii_ffMatrix( &xsize, &ysize, fp ) ;

	if( xsize>1 ){
		fprintf( stderr, "%s: file %s should have a single column \n", __FILE__, filename );
		exit( EXIT_FAILURE );
	}

    gvec=(gsl_vector *)gsl_vector_alloc ( (size_t)ysize );

    for( y=0; y<ysize; y++ ){
        tmp=*(ptr + y );
        gsl_vector_set (gvec, y, tmp);
    }

	fclose(fp);	
	free(ptr);
	return(gvec);
}



/*
===============================
========= main ================
===============================
*/

int main( int argc, char **argv )
{


/*
We need as input:

x  - parameter ensemble matrix
hx - matrix of predicted observations 
R  - observation uncertainty matrix
y  - observation vector

Optional:

hx_bar - vector of predicted observations from mean parameters

Processing:

Calculate Xb from x
Calculate w
Calculate HXb

*/
    
	gsl_matrix *xb ;
	gsl_matrix *hx ;
    gsl_matrix *R ;
    gsl_vector *y ;
	gsl_vector *xa ;
	
	xb =load_ffMatrix_to_gsl_matrix(argv[1]);
	hx =load_ffMatrix_to_gsl_matrix(argv[2]);
	y  =load_ffMatrix_to_gsl_vector(argv[3]);
	R  =load_ffMatrix_to_gsl_matrix(argv[4]);

    xa = fourDEnVar( xb, hx, y, R );

    print_gsl_vector(xa);    

    gsl_matrix_free (xb);
    gsl_vector_free (xa);
    gsl_matrix_free (hx);
    gsl_matrix_free (R);
    gsl_vector_free (y);

	return( EXIT_SUCCESS );
	
}




