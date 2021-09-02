#ifndef FOURDENVAR_ENGINE_H
#define FOURDENVAR_ENGINE_H

#include <gsl_matrix.h>
#include <gsl_blas.h>
#include <gsl_multimin.h>
#include <gsl_linalg.h>

gsl_vector * fourDEnVar( gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_matrix *, gsl_vector * );
gsl_vector * mean_vector_from_matrix( gsl_matrix * );
gsl_matrix * perturbation_matrix( gsl_matrix *, gsl_vector *, float );

void print_gsl_matrix( gsl_matrix * );
void print_gsl_vector( gsl_vector * );

double fourDEnVar_cost_f(const gsl_vector *, void *);
void fourDEnVar_cost_df(const gsl_vector *, void *, gsl_vector *);
void fourDEnVar_cost_fdf (const gsl_vector *, void *, double *, gsl_vector *);

typedef struct  {

    gsl_vector * y; 
    gsl_vector * hx_bar; 
    gsl_matrix * R_inv;
    gsl_matrix * HX_dash_b;

} fourDEnVar_cost_function_vars ;



#endif

