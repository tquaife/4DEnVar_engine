#ifndef FOURDENVAR_ENGINE_H
#define FOURDENVAR_ENGINE_H

#include <gsl_matrix.h>
#include <gsl_blas.h>
#include <gsl_multimin.h>
#include <gsl_linalg.h>

gsl_vector * fourDEnVar( gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_matrix *, gsl_vector * );
gsl_vector * mean_vector_from_matrix( gsl_matrix * );
gsl_matrix * perturbation_matrix( gsl_matrix *, gsl_vector *, float );
gsl_matrix * fourDEnVar_sample_posterior( gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_vector * );

double fourDEnVar_JEval( gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_matrix *, gsl_vector *, gsl_vector * );

void print_gsl_matrix( gsl_matrix * );
void print_gsl_vector( gsl_vector * );

double fourDEnVar_cost_f(const gsl_vector *, void *);
void fourDEnVar_cost_df(const gsl_vector *, void *, gsl_vector *);
void fourDEnVar_cost_fdf (const gsl_vector *, void *, double *, gsl_vector *);

double fourDEnVar_cost_f_error_surface(const gsl_vector *, void *);

typedef struct  {

    gsl_vector * y; 
    gsl_vector * hx_bar; 
    gsl_matrix * R_inv;
    //only need for the new error surfaces:
    gsl_matrix * X_dash_b;
    gsl_matrix * HX_dash_b;

} fourDEnVar_cost_function_vars ;






#endif

