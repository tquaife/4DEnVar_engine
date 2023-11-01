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

x      - parameter ensemble matrix
hx     - matrix of predicted observations 
R      - observation uncertainty matrix
y      - observation vector
hx_bar - vector of predicted observations from mean parameters
x_eval - parameter combinations at which to evaluate to cost function

*/
    
    gsl_matrix *xb ;
    gsl_matrix *hx ;
    gsl_matrix *R ;
    gsl_vector *y ;
    gsl_vector *hx_bar ;
    gsl_matrix *x_eval_matrix ;
    gsl_vector *x_eval_vector ;
    
    double J, tmp;
    int i, j;
    
    xb      =load_ffMatrix_to_gsl_matrix(argv[1]);
    hx      =load_ffMatrix_to_gsl_matrix(argv[2]);
    y       =load_ffMatrix_to_gsl_vector(argv[3]);
    R       =load_ffMatrix_to_gsl_matrix(argv[4]);
    hx_bar  =load_ffMatrix_to_gsl_vector(argv[5]);
   
    x_eval_matrix =load_ffMatrix_to_gsl_matrix(argv[6]);
    x_eval_vector =(gsl_vector *)gsl_vector_alloc(x_eval_matrix->size2);

    for(i=0; i<x_eval_matrix->size1; i++){
        for(j=0; j<x_eval_matrix->size2; j++){
            tmp=gsl_matrix_get(x_eval_matrix,i,j);
            printf("%lf ", tmp);
            gsl_vector_set(x_eval_vector,j,tmp);
        }    
        J = fourDEnVar_JEval( xb, hx, y, R, hx_bar, x_eval_vector );
        printf("%lf\n", J);
    }
    
    gsl_matrix_free (xb);
    gsl_matrix_free (hx);
    gsl_matrix_free (R);
    gsl_vector_free (y);
    gsl_vector_free (hx_bar);
    gsl_matrix_free (x_eval_matrix);    
    gsl_vector_free (x_eval_vector);


    return( EXIT_SUCCESS );
    
}




