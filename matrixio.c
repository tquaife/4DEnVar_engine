#include<matrixio.h>

#define mioVERSION 0.3a

/*

The begginings of a simple library to provide functionality
for reading in vectors and matices, in ascii format. 

*/

/*
** Print message and bomb on error in allocating memory
*/

void	alloc_error( caller_name )
char	*caller_name;
{

	fprintf( stderr, "%s: Error allocating memory\n", caller_name );
	exit( EXIT_FAILURE );
	
}


/*
** ==========================================
** ======== Long int routines =============== 
** ==========================================
*/


/*
** Print a long int matrix to the stdout
*/


void	print_ldMatrix( Matrix_p, x_size, y_size )
long	*Matrix_p, x_size, y_size;
{

	int	i;
	
	for(i=0; i<y_size; ++i )
		print_ldArray( Matrix_p+( x_size*i ), x_size );

	return;


}

/*
** Print a long int array to the stdout
*/

void	print_ldArray( Array_p, size )
long	*Array_p, size;
{

	while( size-- )
		printf( "%ld\t", *Array_p++ );
	
	printf( "\n" );

	return;
}


/*
** Read input file into the matrix array
** Data is strored at the pointer Matrix_p
** and the X and Y coords are stored in
** X_axis and Y_axis (they are calculated 
** by this routine).
*/

long	*read_ascii_ldMatrix( X_size, Y_size )
long	*X_size, *Y_size;
{

	long	element, nlines=0, nelements=0, last_nelements=0;
	long	*Matrix_p;
	
	char	line[ MATRIXIO_MAX_LINE_LEN ];
		
	/* malloc the first byte of memory */
	/* code barfs without this line    */
	if( ( Matrix_p = (long *) malloc( sizeof( long ) ) ) == NULL )
		alloc_error( "read_ascii_ldMatrix" );


	while( fgets( line, MATRIXIO_MAX_LINE_LEN, stdin ) != NULL ){
	
		++nlines;
		nelements=0;
		
		while( get_ldElement( line, &element ) != 0 ){
		
			++nelements;
			
			Matrix_p = ( long *) realloc( ( void *) Matrix_p,\
			      sizeof( long ) * ( ( nlines - 1 ) * last_nelements + nelements ) ) ;
			
			if( Matrix_p == NULL )
				alloc_error( "read_ascii_ldMatrix" );
			 
			*(Matrix_p + ( nlines - 1 ) * last_nelements + nelements -1 ) = element;
		
		}
		
		if( nelements == 0 && nlines > 1 ){
		
			--nlines;
			break;
		
			}else if( nelements == 0 && nlines == 1 ){
			
			--nlines;
		
			}else if( nelements != last_nelements && nlines > 1 ){
			
			fprintf( stderr, "Inconsistant number of columns at line %ld\n", nlines );
			exit( EXIT_FAILURE );
			
			}
		
		
		last_nelements = nelements ;
		
	}

	*X_size = last_nelements ;
	*Y_size = nlines ; 

	return( Matrix_p );

}



/*
As above - but reading from a file.
*/

long	*fread_ascii_ldMatrix( X_size, Y_size, stream )
long	*X_size, *Y_size;
FILE	*stream ;
{

	long	element, nlines=0, nelements=0, last_nelements=0;
	long	*Matrix_p;
	
	char	line[ MATRIXIO_MAX_LINE_LEN ];
		
	/* malloc the first byte of memory */
	/* code barfs without this line    */
	if( ( Matrix_p = (long *) malloc( sizeof( long ) ) ) == NULL )
		alloc_error( "read_ascii_ldMatrix" );


	while( fgets( line, MATRIXIO_MAX_LINE_LEN, stream ) != NULL ){
	
		++nlines;
		nelements=0;
		
		while( get_ldElement( line, &element ) != 0 ){
		
			++nelements;
			
			Matrix_p = ( long *) realloc( ( void *) Matrix_p,\
			      sizeof( long ) * ( ( nlines - 1 ) * last_nelements + nelements ) ) ;
			
			if( Matrix_p == NULL )
				alloc_error( "read_ascii_ldMatrix" );
			 
			*(Matrix_p + ( nlines - 1 ) * last_nelements + nelements -1 ) = element;
		
		}
		
		if( nelements == 0 && nlines > 1 ){
		
			--nlines;
			break;
		
			}else if( nelements == 0 && nlines == 1 ){
			
			--nlines;
		
			}else if( nelements != last_nelements && nlines > 1 ){
			
			fprintf( stderr, "Inconsistant number of columns at line %ld\n", nlines );
			exit( EXIT_FAILURE );
			
			}
		
		
		last_nelements = nelements ;
		
	}

	*X_size = last_nelements ;
	*Y_size = nlines ; 

	return( Matrix_p );

}




/*
** Get elements from an char string, convert them into longs
** and them strip them from the string. Return 0 if none read.
** There's probably a quicker way of doing this using strtol!
*/

char	get_ldElement( Line, Element_p )
char	*Line;
long	*Element_p;
{

	char	first_word[ 20 ];
	int	i=0, j=0;


	for( i=0;i<20;i++ ) first_word[ i ] = ' ';
	i=0 ;


	/* Check that we are not already */
	/* at the end of the line        */

	if( *Line == '\0' ) return( 0 );
	
	/* Skip leading white space and */
	/* return zero if end of string */
	/* is encountered               */
	
	while( isspace( *(Line+i) ) ) 
		if( *(Line+( ++i )) == '\0' ) return( 0 );
	
	
	/* Put the first block of non white space       */
	/* characters into first_word and then          */
	/* convert to a long. Additionally -- If end of */ 
	/* string is encoutered set first character of  */
	/* Line to be the end of string character and   */
	/* return success.                              */
	
	while( ! isspace( *(Line+i) ) ){
		if( ( first_word[j++] = *(Line+(i++) ) ) == '\0' ){
		
			*Element_p = atol( first_word );
			*Line = '\0';
			return( 1 );
		
		}
	}
	
	*Element_p = atol( first_word );
	
	
	/* Copy remainder of Line into the beggining */
	/* of itself (effectively deleting the bit   */
	/* which has just been converted to a long)  */
	
	for( j=0; ( *(Line+j)=*(Line+i) ) != '\0' ; ++j, ++i );
	
	return( 1 );

}





/*
** ==========================================
** ======== Double routines ================= 
** ==========================================
*/


/*
** Print a double matrix to the stdout
*/


void	print_ffMatrix( Matrix_p, x_size, y_size )
long	x_size, y_size;
double	*Matrix_p;
{

	int	i;
	
	for(i=0; i<y_size; ++i )
		print_ffArray( Matrix_p+( x_size*i ), x_size );

	return;


}

/*
** Print a double precision array to the stdout
*/

void	print_ffArray( Array_p, size )
long	size;
double	*Array_p; 
{

	while( size-- )
		printf( "%f\t", *Array_p++ );
	
	printf( "\n" );

	return;
}


/*
** Double version: 
** Read input file (from stdin) into the 
** matrix array.
** Data is strored at the pointer Matrix_p
** and the X and Y coords are stored in
** X_axis and Y_axis (they are calculated 
** by this routine).
*/

double	*read_ascii_ffMatrix( X_size, Y_size )
long	*X_size, *Y_size;
{

	long	nlines=0, nelements=0, last_nelements=0;
	double	*Matrix_p, element;
	
	char	line[ MATRIXIO_MAX_LINE_LEN ];
		
	/* malloc the first byte of memory */
	/* code barfs without this line    */
	if( ( Matrix_p = ( double *) malloc( sizeof( double ) ) ) == NULL )
		alloc_error( "read_ascii_ffMatrix" );


	while( fgets( line, MATRIXIO_MAX_LINE_LEN, stdin ) != NULL ){
	
		++nlines;
		nelements=0;
		
		while( get_ffElement( line, &element ) != 0 ){
		
			++nelements;
			
			Matrix_p = ( double *) realloc( ( void *) Matrix_p,\
			      sizeof( double ) * ( ( nlines - 1 ) * last_nelements + nelements ) ) ;
			
			if( Matrix_p == NULL )
				alloc_error( "read_ascii_ffMatrix" );
			 
			*(Matrix_p + ( nlines - 1 ) * last_nelements + nelements -1 ) = element;
		
		}
		
		if( nelements == 0 && nlines > 1 ){
		
			--nlines;
			break;
		
			}else if( nelements == 0 && nlines == 1 ){
			
			--nlines;
		
			}else if( nelements != last_nelements && nlines > 1 ){
			
			fprintf( stderr, "Inconsistant number of columns at line %ld\n", nlines );
			exit( EXIT_FAILURE );
			
			}
		
		
		last_nelements = nelements ;
		
	}

	*X_size = last_nelements ;
	*Y_size = nlines ; 

	return( Matrix_p );

}


/* Double version:
** Get elements from an char string, convert them into longs
** and them strip them from the string. Return 0 if none read.
** There's probably a quicker way of doing this using strtol!
*/

char	get_ffElement( Line, Element_p )
char	*Line;
double	*Element_p;
{

	char	first_word[ 20 ];
	int	i=0, j=0;


	for( i=0;i<20;i++ ) first_word[ i ] = ' ';
	i=0 ;

	/* Check that we are not already */
	/* at the end of the line        */

	if( *Line == '\0' ) return( 0 );
	
	/* Skip leading white space and */
	/* return zero if end of string */
	/* is encountered               */
	
	while( isspace( *(Line+i) ) ) 
		if( *(Line+( ++i )) == '\0' ) return( 0 );
	
	
	/* Put the first block of non white space       */
	/* characters into first_word and then          */
	/* convert to a long. Additionally -- If end of */ 
	/* string is encoutered set first character of  */
	/* Line to be the end of string character and   */
	/* return success.                              */
	
	while( ! isspace( *(Line+i) ) ){
		if( ( first_word[j++] = *(Line+(i++) ) ) == '\0' ){
		
			*Element_p = atof( first_word );
			*Line = '\0';
			return( 1 );
		
		}
	}
	
	*Element_p = atof( first_word );
	
	
	/* Copy remainder of Line into the beggining */
	/* of itself (effectively deleting the bit   */
	/* which has just been converted to a long)  */
	
	for( j=0; ( *(Line+j)=*(Line+i) ) != '\0' ; ++j, ++i );
	
	return( 1 );

}



/*
** Double version:
** Read input file into the matrix array
** Data is strored at the pointer Matrix_p
** and the X and Y coords are stored in
** X_axis and Y_axis (they are calculated 
** by this routine).
*/

double	*fread_ascii_ffMatrix( X_size, Y_size, stream )
long	*X_size, *Y_size;
FILE	*stream;
{

	long	nlines=0, nelements=0, last_nelements=0;
	double	*Matrix_p, element;
	
	char	line[ MATRIXIO_MAX_LINE_LEN ];
		
	/* malloc the first byte of memory */
	/* code barfs without this line    */
	if( ( Matrix_p = ( double *) malloc( sizeof( double ) ) ) == NULL )
		alloc_error( "fread_ascii_ffMatrix" );


	while( fgets( line, MATRIXIO_MAX_LINE_LEN, stream ) != NULL ){
	
		++nlines;
		nelements=0;
		
		while( get_ffElement( line, &element ) != 0 ){
		
			++nelements;
			
			Matrix_p = ( double *) realloc( ( void *) Matrix_p,\
			      sizeof( double ) * ( ( nlines - 1 ) * last_nelements + nelements ) ) ;
			
			if( Matrix_p == NULL )
				alloc_error( "fread_ascii_ffMatrix" );
			 
			*(Matrix_p + ( nlines - 1 ) * last_nelements + nelements -1 ) = element;
		
		}
		
		if( nelements == 0 && nlines > 1 ){
		
			--nlines;
			break;
		
			}else if( nelements == 0 && nlines == 1 ){
			
			--nlines;
		
			}else if( nelements != last_nelements && nlines > 1 ){
			
			fprintf( stderr, "Inconsistant number of columns at line %ld\n", nlines );
			exit( EXIT_FAILURE );
			
			}
		
		
		last_nelements = nelements ;
		
	}

	*X_size = last_nelements ;
	*Y_size = nlines ; 

	return( Matrix_p );

}




/*
** ==========================================
** ========== Float routines =============== 
** ==========================================
*/


/*
** Print a single precision matrix to the stdout
*/


void	print_fMatrix( Matrix_p, x_size, y_size )
long	x_size, y_size;
float	*Matrix_p;
{

	int	i;
	
	for(i=0; i<y_size; ++i )
		print_fArray( Matrix_p+( x_size*i ), x_size );

	return;


}

/*
** Print a single precision array to the stdout
*/

void	print_fArray( Array_p, size )
long	size;
float	*Array_p; 
{

	while( size-- )
		printf( "%f\t", *Array_p++ );
	
	printf( "\n" );

	return;
}


/*
** single precision version: 
** Read input file (from stdin) into the 
** matrix array.
** Data is strored at the pointer Matrix_p
** and the X and Y coords are stored in
** X_axis and Y_axis (they are calculated 
** by this routine).
*/

float	*read_ascii_fMatrix( X_size, Y_size )
long	*X_size, *Y_size;
{

	long	nlines=0, nelements=0, last_nelements=0;
	float	*Matrix_p, element;
	
	char	line[ MATRIXIO_MAX_LINE_LEN ];
		
	/* malloc the first byte of memory */
	/* code barfs without this line    */
	if( ( Matrix_p = ( float *) malloc( sizeof( float ) ) ) == NULL )
		alloc_error( "read_ascii_fMatrix" );


	while( fgets( line, MATRIXIO_MAX_LINE_LEN, stdin ) != NULL ){
	
		++nlines;
		nelements=0;
		
		while( get_fElement( line, &element ) != 0 ){
		
			++nelements;
			
			Matrix_p = ( float *) realloc( ( void *) Matrix_p,\
			      sizeof( float ) * ( ( nlines - 1 ) * last_nelements + nelements ) ) ;
			
			if( Matrix_p == NULL )
				alloc_error( "read_ascii_fMatrix" );
			 
			*(Matrix_p + ( nlines - 1 ) * last_nelements + nelements -1 ) = element;
		
		}
		
		if( nelements == 0 && nlines > 1 ){
		
			--nlines;
			break;
		
			}else if( nelements == 0 && nlines == 1 ){
			
			--nlines;
		
			}else if( nelements != last_nelements && nlines > 1 ){
			
			fprintf( stderr, "Inconsistant number of columns at line %ld\n", nlines );
			exit( EXIT_FAILURE );
			
			}
		
		
		last_nelements = nelements ;
		
	}

	*X_size = last_nelements ;
	*Y_size = nlines ; 

	return( Matrix_p );

}


/* Float version:
** Get elements from an char string, convert them into longs
** and them strip them from the string. Return 0 if none read.
** There's probably a quicker way of doing this using strtol!
*/

char	get_fElement( Line, Element_p )
char	*Line;
float	*Element_p;
{

	char	first_word[ 20 ];
	int	i=0, j=0;
	
	
	for( i=0;i<20;i++ ) first_word[ i ] = ' ';
	i=0 ;


	/* Check that we are not already */
	/* at the end of the line        */

	if( *Line == '\0' ) return( 0 );
	
	/* Skip leading white space and */
	/* return zero if end of string */
	/* is encountered               */
	
	while( isspace( *(Line+i) ) ) 
		if( *(Line+( ++i )) == '\0' ) return( 0 );
	
	
	/* Put the first block of non white space       */
	/* characters into first_word and then          */
	/* convert to a long. Additionally -- If end of */ 
	/* string is encoutered set first character of  */
	/* Line to be the end of string character and   */
	/* return success.                              */
	
	while( ! isspace( *(Line+i) ) ){
		if( ( first_word[j++] = *(Line+(i++) ) ) == '\0' ){
		
			*Element_p = (float) atof( first_word );
			*Line = '\0';
			return( 1 );
		
		}
	}
	
	*Element_p = (float) atof( first_word );
	
	
	/* Copy remainder of Line into the beggining */
	/* of itself (effectively deleting the bit   */
	/* which has just been converted to a long)  */
	
	for( j=0; ( *(Line+j)=*(Line+i) ) != '\0' ; ++j, ++i );
	
	return( 1 );

}



/*
** Float version:
** Read input file into the matrix array
** Data is strored at the pointer Matrix_p
** and the X and Y coords are stored in
** X_axis and Y_axis (they are calculated 
** by this routine).
*/

float	*fread_ascii_fMatrix( X_size, Y_size, stream )
long	*X_size, *Y_size;
FILE	*stream;
{

	long	nlines=0, nelements=0, last_nelements=0;
	float	*Matrix_p, element;
	
	char	line[ MATRIXIO_MAX_LINE_LEN ];
		
	/* malloc the first byte of memory */
	/* code barfs without this line    */
	if( ( Matrix_p = ( float *) malloc( sizeof( float ) ) ) == NULL )
		alloc_error( "fread_ascii_fMatrix" );


	while( fgets( line, MATRIXIO_MAX_LINE_LEN, stream ) != NULL ){
	
		++nlines;
		nelements=0;
		
		while( get_fElement( line, &element ) != 0 ){
		
			++nelements;
			
			Matrix_p = ( float *) realloc( ( void *) Matrix_p,\
			      sizeof( float ) * ( ( nlines - 1 ) * last_nelements + nelements ) ) ;
			
			if( Matrix_p == NULL )
				alloc_error( "sread_ascii_fMatrix" );
			 
			*(Matrix_p + ( nlines - 1 ) * last_nelements + nelements -1 ) = element;
		
		}
		
		if( nelements == 0 && nlines > 1 ){
		
			--nlines;
			break;
		
			}else if( nelements == 0 && nlines == 1 ){
			
			--nlines;
		
			}else if( nelements != last_nelements && nlines > 1 ){
			
			fprintf( stderr, "Inconsistant number of columns at line %ld\n", nlines );
			exit( EXIT_FAILURE );
			
			}
		
		
		last_nelements = nelements ;
		
	}

	*X_size = last_nelements ;
	*Y_size = nlines ; 

	return( Matrix_p );

}





/*
** Flip a floating point matrix to exchange
** rows and columns. This is needed as the
** matrix format in this library is "flipped"
** compared to the LAPACK libraries. 
*/

void	flip_fMatrix( Matrix_p, x_size, y_size )
float	*Matrix_p;
long	y_size, x_size;
{

	long	i=0, j=0, k=0;
	float	*tmp_matrix;
	
	/*allocate memory to temporary matrix*/
	
	if( ( tmp_matrix = ( float *) malloc( x_size * y_size * sizeof( float ) ) ) == NULL )
		alloc_error( "flip_fMatrix" );

	
	/*copy matrix into tmp*/
	
	for( i=0; i<x_size*y_size; ++i )
		*( tmp_matrix + i ) = *( Matrix_p + i );
	
	
	for( i=0; i<x_size; i++ ){
		for( j=0; j<y_size; j++ ){
			
			*( Matrix_p + k++ )= *( tmp_matrix + i + ( j * x_size ) ); 
		
			}
		}
		
	free( tmp_matrix );

	return;
}





