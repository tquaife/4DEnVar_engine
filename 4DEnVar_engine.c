#include<4DEnVar_engine.h>
#include<math.h>

gsl_vector * fourDEnVar( gsl_matrix * xb, gsl_matrix * hx, gsl_vector * y, gsl_matrix * R, gsl_vector * hx_bar )
/*
An implementation of 4DEnVar
*/
{

    gsl_vector *xb_bar ;
    //gsl_vector *hx_bar ;
    gsl_matrix *X_dash_b ;
    gsl_matrix *HX_dash_b ;
    gsl_vector *w = gsl_vector_calloc(xb->size2) ;  /*calloc ensures w=0*/
    int nens=xb->size1;
    int status ;
    size_t iter=0;
    float scale ;
        
    fourDEnVar_cost_function_vars cost_vars;
    
    const gsl_multimin_fdfminimizer_type *solver_type;
    gsl_multimin_fdfminimizer *solver;
    gsl_multimin_function_fdf cost_func;

    int signum;
    gsl_permutation *p=gsl_permutation_alloc(R->size1);
    gsl_matrix *R_inv = gsl_matrix_calloc(R->size1, R->size2);

    solver_type = gsl_multimin_fdfminimizer_vector_bfgs2;
    solver = gsl_multimin_fdfminimizer_alloc (solver_type, w->size);

    /*invert the R matrix*/
    gsl_linalg_LU_decomp(R, p, &signum);
    gsl_linalg_LU_invert(R, p, R_inv);
        
    /*calculate the mean of each parameter 
    in the ensemble*/
    xb_bar = mean_vector_from_matrix(xb);
    //print_gsl_vector(xb_bar);
    //printf("***\n");


    /*mean of modelled observations*/
    /*NOTE - should be supplied separately as h(x)*/
    //hx_bar = mean_vector_from_matrix(hx);
    //print_gsl_vector(hx_bar);
    //printf("***\n");

    
    /*calculate the perturbation matrix
    eqn 21 in Pinnington 2020*/
    scale=1./sqrt(nens-1);
    X_dash_b = perturbation_matrix(xb,xb_bar,scale);

    /*calculate HXb matrix
    eqn 26 in Pinnington 2020
    */
    HX_dash_b = perturbation_matrix(hx,hx_bar,scale);
    //print_gsl_matrix(HX_dash_b);


    /*n.b. as we are setting x0=xb_bar then w=0 so
    we do not need to compute the initial transformation
    into ensemble space as w is initialised as 0
        
    As a consquence, we don't need an x0 variable. Obvs.
    */


    /*set up structure holding everything
    needed to calculate the cost function #
    (other than w)*/
    cost_vars.R_inv=R_inv;
    cost_vars.HX_dash_b=HX_dash_b;
    cost_vars.hx_bar=hx_bar;
    cost_vars.y=y;

    /*set up the cost function structure*/
    cost_func.n = w->size;
    cost_func.f = fourDEnVar_cost_f;
    cost_func.df = fourDEnVar_cost_df;
    cost_func.fdf = fourDEnVar_cost_fdf;
    cost_func.params = &cost_vars;

    /*initialise the minimiser*/
    gsl_multimin_fdfminimizer_set(solver, &cost_func, w, 0.01, 1e-10);

    /*iterate the minimiser*/
    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(solver);
        if (status) break;
        status = gsl_multimin_test_gradient(solver->gradient, 1e-16);
    } while (status == GSL_CONTINUE && iter < 100);

    //printf("## n_iter=%ld stat=%d (enoprog=%d)\n", iter, status, GSL_ENOPROG);

    /*translate w back into observation space*/
    /*n.b. writing the answer (xa) over xb_bar as we need to 
    add that back in anyway, so convenient given dgemv */
    gsl_blas_dgemv(CblasNoTrans, 1.0,X_dash_b, solver->x, 1.0, xb_bar);

    return(xb_bar);
}


double fourDEnVar_cost_f(const gsl_vector *w, void *p)
{

    double bgrnd_term ; 
    double obs_term ; 

    double J;
    
    const gsl_vector *hx_bar    = ((fourDEnVar_cost_function_vars *) p)->hx_bar;
    const gsl_vector *y         = ((fourDEnVar_cost_function_vars *) p)->y;
    const gsl_matrix *R_inv     = ((fourDEnVar_cost_function_vars *) p)->R_inv;
    const gsl_matrix *HX_dash_b = ((fourDEnVar_cost_function_vars *) p)->HX_dash_b;
    
    gsl_vector * work1=gsl_vector_alloc(y->size) ;
    gsl_vector * work2=gsl_vector_alloc(y->size) ;

    /* w^Tw */
    gsl_blas_ddot(w, w, &bgrnd_term);

    /* hx-y */    
    gsl_vector_memcpy(work1,hx_bar);
    gsl_blas_daxpy(-1.,y,work1);    

    /* Hx'b * w + (hx-y)*/   
    /* work1 should be hx-y*/
    gsl_blas_dgemv(CblasNoTrans, 1.0, HX_dash_b, w, 1.0, work1);

    /* R^-1*(HX'b*w+hx-y) */
    /* work1 should now be (HX'b*w+hx-y)*/
    gsl_vector_memcpy(work2,work1);
    gsl_blas_dgemv(CblasTrans, 1.0, R_inv, work2, 0.0, work1);
    
    /*full observation error term*/
    /* work1 should now be R^-1*(HX'b*w+hx-y)*/
    /* work2 should now be (HX'b*w+hx-y)*/
    gsl_blas_ddot(work1, work2, &obs_term);
    
    J=0.5*bgrnd_term+0.5*obs_term;    
    return(J);    
}

void fourDEnVar_cost_df(const gsl_vector *w, void *p, gsl_vector *df)
{

    const gsl_vector *hx_bar    = ((fourDEnVar_cost_function_vars *) p)->hx_bar;
    const gsl_vector *y         = ((fourDEnVar_cost_function_vars *) p)->y;
    const gsl_matrix *R_inv     = ((fourDEnVar_cost_function_vars *) p)->R_inv;
    const gsl_matrix *HX_dash_b = ((fourDEnVar_cost_function_vars *) p)->HX_dash_b;
    
    gsl_vector * work1=gsl_vector_alloc(y->size) ;
    gsl_vector * work2=gsl_vector_alloc(y->size) ;


    /* hx-y */    
    gsl_vector_memcpy(work1,hx_bar);
    gsl_blas_daxpy(-1.,y,work1);    

    /*Hx'b * w + (hx-y)*/   
    /* work1 should be hx-y*/
    gsl_blas_dgemv(CblasNoTrans, 1.0, HX_dash_b, w, 1.0, work1);

    /* R^-1*(HX'b*w+hx-y) */
    /* work1 should now be (HX'b*w+hx-y)*/
    gsl_vector_memcpy(work2,work1);
    gsl_blas_dgemv(CblasTrans, 1.0, R_inv, work2, 0.0, work1);
    
    /*full gradient vector for cost function*/
    /* work1 should now be R^-1*(HX'b*w+hx-y)*/
    /*copy w into df and add it onto the matrix-vector
    product (HX'b)^T * R^-1*(HX'b*w+hx-y) */
    gsl_vector_memcpy(df,w);
    gsl_blas_dgemv(CblasTrans, 1.0, HX_dash_b, work1, 1.0, df);
   
    /*df should now be: 
        w + (HX'b)^T*R^-1*(HX'b*w+hx-y) 
      so we are done...*/
       
    return;
}


void fourDEnVar_cost_fdf (const gsl_vector *w, void *p, double *f, gsl_vector *df)
{
  *f = fourDEnVar_cost_f(w, p);
  fourDEnVar_cost_df(w, p, df);
  return;
}


gsl_vector * mean_vector_from_matrix( gsl_matrix * gmat )
/*
Compute means across rows of a matrix and return as a GSL vector 
*/
{

    int x, y;
    float tmp ;
    gsl_vector *gvec = gsl_vector_alloc(gmat->size1);

    for( y=0; y<gmat->size1; y++ ){    
        tmp=0.0;
        for( x=0; x<gmat->size2; x++ ){
            tmp+=gsl_matrix_get(gmat, y, x );
        }
        gsl_vector_set( gvec, y, tmp/(float)gmat->size2 );
    }
    return(gvec);

}



gsl_matrix * perturbation_matrix( gsl_matrix * gmat, gsl_vector * gvec, float scale )
/*
Calculate a perturbation matrix by subtract the ith element of
gvec from each element of the ith row of gmat.
gvec should be the mean of each row of gmat, for example
calculated by mean_vector_from_matrix (but that is not done
inside this function to allow some flexibility)
*/
{
    int x,y;
    float tmp;
    gsl_matrix *pert = gsl_matrix_alloc(gmat->size1,gmat->size2);
    for( x=0; x<gmat->size2; x++ ){
        for( y=0; y<gmat->size1; y++ ){
            tmp=gsl_matrix_get(gmat,y,x)-gsl_vector_get(gvec,y);
            gsl_matrix_set (pert, y, x, tmp*scale);
        }
    }
    return(pert);
}


void print_gsl_matrix( gsl_matrix * gmat )
/*
Print an ascii representation of the GSL matrix gmat to the stdout
*/
{
    int x, y ;
    double tmp;

    for( y=0; y<gmat->size1; y++ ){    
        for( x=0; x<gmat->size2; x++ ){
            tmp=gsl_matrix_get(gmat, y, x );
            printf("%lf ",tmp); 
        }
        printf("\n");
    }
    return;
}


void print_gsl_vector( gsl_vector * gvec )
/*
Print an ascii representation of the GSL matrix gmat to the stdout
*/
{
    int y ;
    double tmp;

    for( y=0; y<gvec->size; y++ ){
        tmp=gsl_vector_get(gvec, y);
        printf("%lf\n",tmp); 
    }
    return;
}

