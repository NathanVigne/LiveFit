/*
#include <gsl/gsl_statistics.h>
#include <stdio.h>

int main(void)
{
    double data[5] = {17.2, 18.1, 16.5, 18.3, 12.6};
    double mean, variance, largest, smallest;

    mean = gsl_stats_mean(data, 1, 5);
    variance = gsl_stats_variance(data, 1, 5);
    largest = gsl_stats_max(data, 1, 5);
    smallest = gsl_stats_min(data, 1, 5);

    printf("The dataset is %g, %g, %g, %g, %g\n", data[0], data[1], data[2], data[3], data[4]);

    printf("The sample mean is %g\n", mean);
    printf("The estimated variance is %g\n", variance);
    printf("The largest value is %g\n", largest);
    printf("The smallest value is %g\n", smallest);
    return 0;
}
*/
/*
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <stdlib.h>


#define N 100      // number of data points to fit 

#define TMAX (3.0) // time variable in [0,TMAX] 

struct data
{
    size_t n;
    double *t;
    double *y;
};

int expb_f(const gsl_vector *x, void *data, gsl_vector *f)
{
    size_t n = ((struct data *) data)->n;
    double *t = ((struct data *) data)->t;
    double *y = ((struct data *) data)->y;

    double A = gsl_vector_get(x, 0);
    double lambda = gsl_vector_get(x, 1);
    double b = gsl_vector_get(x, 2);

    size_t i;

    for (i = 0; i < n; i++) {
        // Model Yi = A * exp(-lambda * t_i) + b 
        double Yi = A * exp(-lambda * t[i]) + b;
        gsl_vector_set(f, i, Yi - y[i]);
    }

    return GSL_SUCCESS;
}

int expb_df(const gsl_vector *x, void *data, gsl_matrix *J)
{
    size_t n = ((struct data *) data)->n;
    double *t = ((struct data *) data)->t;

    double A = gsl_vector_get(x, 0);
    double lambda = gsl_vector_get(x, 1);

    size_t i;

    for (i = 0; i < n; i++) {
        // Jacobian matrix J(i,j) = dfi / dxj, 
        // where fi = (Yi - yi)/sigma[i],      
        //       Yi = A * exp(-lambda * t_i) + b
        // and the xj are the parameters (A,lambda,b) 
        double e = exp(-lambda * t[i]);
        gsl_matrix_set(J, i, 0, e);
        gsl_matrix_set(J, i, 1, -t[i] * A * e);
        gsl_matrix_set(J, i, 2, 1.0);
    }

    return GSL_SUCCESS;
}

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);
    double rcond;

    // compute reciprocal condition number of J(x) 
    gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stderr,
            "iter %2zu: A = %.4f, lambda = %.4f, b = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
            iter,
            gsl_vector_get(x, 0),
            gsl_vector_get(x, 1),
            gsl_vector_get(x, 2),
            1.0 / rcond,
            gsl_blas_dnrm2(f));
}

int main(void)
{
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    const size_t n = N;
    const size_t p = 3;

    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc(p, p);
    double t[N], y[N], weights[N];
    struct data d = {n, t, y};
    double x_init[3] = {1.0, 1.0, 0.0}; // starting values 
    gsl_vector_view x = gsl_vector_view_array(x_init, p);
    gsl_vector_view wts = gsl_vector_view_array(weights, n);
    gsl_rng *r;
    double chisq, chisq0;
    int status, info;
    size_t i;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;

    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    // define the function to be minimized 
    fdf.f = expb_f;
    fdf.df = expb_df;  //set to NULL for finite-difference Jacobian 
    fdf.fvv = NULL;   // not using geodesic acceleration 
    fdf.n = n;
    fdf.p = p;

    fdf.params = &d;

    // this is the data to be fitted 
    for (i = 0; i < n; i++) {
        double ti = i * TMAX / (n - 1.0);
        double yi = 1.0 + 5 * exp(-1.5 * ti);
        double si = 0.1 * yi;
        double dy = gsl_ran_gaussian(r, si);

        t[i] = ti;
        y[i] = yi + dy;
        weights[i] = 1.0 / (si * si);
        printf("data: %g %g %g\n", ti, y[i], si);
    };

    // allocate workspace with default parameters 
    w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);

    // initialize solver with starting point and weights 
    gsl_multifit_nlinear_winit(&x.vector, &wts.vector, &fdf, w);

    // compute initial cost function 
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    // solve the system with a maximum of 100 iterations 
    status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol, callback, NULL, &info, w);

    // compute covariance of best fit parameters 
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar(J, 0.0, covar);

    // compute final cost 
    gsl_blas_ddot(f, f, &chisq);

#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar, i, i))

    fprintf(stderr,
            "summary from method '%s/%s'\n",
            gsl_multifit_nlinear_name(w),
            gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));

    {
        double dof = n - p;
        double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

        fprintf(stderr, "A      = %.5f +/- %.5f\n", FIT(0), c * ERR(0));
        fprintf(stderr, "lambda = %.5f +/- %.5f\n", FIT(1), c * ERR(1));
        fprintf(stderr, "b      = %.5f +/- %.5f\n", FIT(2), c * ERR(2));
    }

    fprintf(stderr, "status = %s\n", gsl_strerror(status));

    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);
    gsl_rng_free(r);

    return 0;}*/
/*
#include <gsl/gsl_linalg.h>
#include <iostream>

int main()
{
    double A_data[] = {0.57092943,
                       0.00313503,
                       0.88069151,
                       0.39626474,
                       0.33336008,
                       0.01876333,
                       0.12228647,
                       0.40085702,
                       0.55534451,
                       0.54090141,
                       0.85848041,
                       0.62154911,
                       0.64111484,
                       0.8892682,
                       0.58922332,
                       0.32858322};

    double b_data[] = {1.5426693, 0.74961678, 2.21431998, 2.14989419};

    // Access the above C arrays through GSL views
    gsl_matrix_view A = gsl_matrix_view_array(A_data, 4, 4);
    gsl_vector_view b = gsl_vector_view_array(b_data, 4);

    // Print the values of A and b using GSL print functions
    std::cout << "A = \n";
    gsl_matrix_fprintf(stdout, &A.matrix, "%lf");

    std::cout << "\nb = \n";
    gsl_vector_fprintf(stdout, &b.vector, "%lf");

    // Allocate memory for the solution vector x and the permutation perm:
    gsl_vector *x = gsl_vector_alloc(4);
    gsl_permutation *perm = gsl_permutation_alloc(4);

    // Decompose A into the LU form:
    int signum;
    gsl_linalg_LU_decomp(&A.matrix, perm, &signum);

    // Solve the linear system
    gsl_linalg_LU_solve(&A.matrix, perm, &b.vector, x);

    // Print the solution
    std::cout << "\nx = \n";
    gsl_vector_fprintf(stdout, x, "%lf");

    // Release the memory previously allocated for x and perm
    gsl_vector_free(x);
    gsl_permutation_free(perm);}*/

#include "mainwidget.h"

#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>

QT_USE_NAMESPACE

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWidget w;
    w.resize(720, 480);
    w.show();

    return a.exec();
}
