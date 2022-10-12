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
    int nens=xb->size2;
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
    scale=1./sqrt((float)nens-1.);
    X_dash_b = perturbation_matrix(xb,xb_bar,scale);

    /*calculate HXb matrix
    eqn 26 in Pinnington 2020
    */
    HX_dash_b = perturbation_matrix(hx,hx_bar,scale);
    //print_gsl_matrix(HX_dash_b);


    /*n.b. as we are setting x0=xb_bar then w=0 so
    we do not need to compute the initial transformation
    into ensemble space as w is initialised as 0.        
    As a consequence, we don't need an x0 variable. Obvs.
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


gsl_matrix * fourDEnVar_sample_posterior( gsl_matrix * xb, gsl_matrix * hx, gsl_matrix * R, gsl_vector * hx_bar, gsl_vector *xa )
{

    gsl_vector *xb_bar ;
    gsl_matrix *X_dash_b ;
    gsl_matrix *HX_dash_b ;
    int nens=xb->size2, i=0, j=0;
    float scale,tmp;

    int signum_r;
    gsl_permutation *p_r=gsl_permutation_alloc(R->size1);
    gsl_matrix *R_inv = gsl_matrix_calloc(R->size1, R->size2);
    int signum_w;
    gsl_permutation *p_w=gsl_permutation_alloc(hx->size2);
    gsl_matrix *Wa = gsl_matrix_calloc(hx->size2, hx->size2);

    gsl_matrix *work1 = gsl_matrix_alloc(hx->size1, hx->size2);
    gsl_matrix *work2 = gsl_matrix_calloc(hx->size2, hx->size2);

    gsl_matrix *X_dash_a = gsl_matrix_alloc(xb->size1, xb->size2);
    gsl_matrix *X_a = gsl_matrix_alloc(xb->size1, xb->size2);


    /*set I*/
    gsl_matrix_set_identity(work2);

    /*invert the R matrix*/
    gsl_linalg_LU_decomp(R, p_r, &signum_r);
    gsl_linalg_LU_invert(R, p_r, R_inv);
    
    /*calculate the mean of each parameter 
    in the ensemble*/
    xb_bar = mean_vector_from_matrix(xb);
    
    /*calculate the perturbation matrix
    eqn 21 in Pinnington 2020
    
    Note - that paper applies a scaling of:
    scale=1./sqrt((float)nens-1.);
    but that is error. Here using scale=1.
    See also step (2) below.
    */
    scale=1.;
    X_dash_b = perturbation_matrix(xb,xb_bar,scale);

    /*calculate HXb matrix
    eqn 26 in Pinnington 2020*/
    HX_dash_b = perturbation_matrix(hx,hx_bar,scale);

    /*Wa=sqrt(I+Y'b^T*R^-1*Y'b)
    See eqn A16 in Pinnington et al. 2021
    Note that Y'd == HX'b due to a change 
    in nomenclature between the two papers
    */

    /*1. work1=R^-1*HX'b*/
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, R_inv, HX_dash_b, 0.0, work1);

    /*2. (HX'b)^T*work1 + I
    Note: work2 = I
    Note: use use of scale here differs from Pinnington et al. 2021
    the equations in the appendix of that paper appear to be incorrect
    */
    
    scale=1./((float)nens-1.); 
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, scale, HX_dash_b, work1, 1.0, work2);

    /* 3. work2 now contains I+Y'b^T*R^-1*Y'b
    so need to take square root to give Wa. 
    Need to zero the upper left triangle as GSL
    leaves the original matrix in there.
    
    (n.b. tested the zeroing of the UL by
    confirming it was necessary to make A=LL^T)
    */

    gsl_linalg_cholesky_decomp1(work2);  
    for(i=0;i<work2->size1;i++)
        for(j=i+1;j<work2->size2;j++)
           gsl_matrix_set(work2, i, j, 0.0);
            
    /*3a. Invert the square root matrix*/
    gsl_linalg_LU_decomp(work2, p_w, &signum_w);
    gsl_linalg_LU_invert(work2, p_w, Wa);

        
    /*4. X'a = X'b * Wa
    Note - work2 contains Wa.
    */     

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X_dash_b, Wa, 0.0, X_dash_a);
    
    /*Now compute Xa=X'a+xa
    */
    for( i=0; i<X_dash_a->size2; i++ ){
        for( j=0; j<X_dash_a->size1; j++ ){
            tmp=gsl_matrix_get(X_dash_a,j,i)+gsl_vector_get(xa,j);
            gsl_matrix_set (X_a, j, i, tmp);
        }
    }
    
    return(X_a);
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

