
#include <iostream>
#include <unistd.h>
#include <cstring>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

const double ALPHA = 0.01;
const double BETA  = 10;
const int    DIMENSION = 10;
const double EPS   = 1E-4;

typedef double* (*matvec_t)(double*, const double* const);

double*         init_matrix();
gsl_vector*     eigenvalues(double *matrix);
inline double   eigenvalue_max(gsl_vector* eigvals);
inline double   eigenvalue_min(gsl_vector* eigvals);

double* copyvec(cosnt double* const vec);
double* normvec(const double* const vec);
double* subvec(double* res, const double* const from, const double* const what);
double* matvec(double* __restrict res, const double* const __restrict vec);
double  SI_solver(double* __restrict f, const double* __restrict x0, matvec_t matvec);

// =========================================================================

int main()
{
    double* matrix    = init_matrix();
    gsl_vector* vec   = eigenvalues(matrix);
    double eigval_max = eigenvalue_max(vec);
    double eigval_min = eigenvalue_min(vec);

    // for(int i = 0; i < DIMENSION; i++)
    //     printf("%lf\n", gsl_vector_get(vec, i));

    printf("min = %lg, max = %lg\n", eigval_max, eigval_min);

    SI_solver(NULL, NULL, matvec);

    gsl_vector_free(vec);
}

// =========================================================================

double*     init_matrix()
{
    double* matrix;
    try
    {
        matrix = new double [DIMENSION*DIMENSION];
    }
    catch (const std::bad_alloc& ex)
    {
        return NULL;
    }

    matrix[0] = 2 + BETA;
    for(int i = 1; i < DIMENSION; i++)
        matrix[i*DIMENSION + i] = 2 + ALPHA;
    for(int i = 1; i < DIMENSION - 1; i++)
        matrix[i*DIMENSION + i + 1] = -1;
    for(int i = 1; i < DIMENSION; i++)
        matrix[i*DIMENSION + i - 1] = -1;

    return matrix;
}

// Returns vector sorted with eigenvalues already sorted
// WARNING:vector needs to be deleted manually
gsl_vector* eigenvalues(double *mat)
{
    gsl_matrix_view m = gsl_matrix_view_array(mat, DIMENSION,
                                                   DIMENSION);
    gsl_vector *eval  = gsl_vector_alloc(DIMENSION);
    gsl_matrix *evec  = gsl_matrix_alloc(DIMENSION, DIMENSION);

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(DIMENSION);
    gsl_eigen_symmv(&m.matrix, eval, evec, w);
    gsl_eigen_symmv_free(w);

    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    gsl_matrix_free(evec);
    return eval;
}

// Expects vector already sorted in abs descending order
inline double eigenvalue_max(gsl_vector* eigvals)
{
    double eigeval = gsl_vector_get(eigvals, 0);
    return std::max(eigeval, -eigeval);
}

// Expects vector already sorted in abs descending order
inline double eigenvalue_min(gsl_vector* eigvals)
{
    double eigeval = gsl_vector_get(eigvals, DIMENSION - 1);
    return std::max(eigeval, -eigeval);
}

// Allocates memory, free it by yourself.
// Returns NULL in case new fails.
double* copyvec(cosnt double* const vec)
{
    double* x;
    try
    {
        x = new double [DIMENSION];
    }
    catch(const std::bad_alloc& ex)
    {
        return NULL;
    }

    return memcpy(x, vec, DIMENSION);
}

double* normvec(const double* const vec)
{
    double sum = 0;
    for(int i = 0; i < DIMENSION; i++)
        sum += vec[i]*vec[i];
    return sqrt(sum);
}

// Subtracts 'what' from 'from' and stores result in res.
double* subvec(double* res = NULL, const double* const from, const double* const what)
{
    if (res == NULL)    res = from;

    for (int i = 0; i < DIMENSION; i++)
        res[i] = from[i] - what[i];
    return res;
}

// Function only works with given matrix
double* matvec(double* __restrict res, const double* const __restrict vec)
{
    res[0] = (2 + BETA)*vec[0] - vec[1];

    for(int i = 1; i < DIMENSION - 1; i++)
        res[i] = -vec[i-1] + (2 + ALPHA)*vec[i] - vec[i+1];

    res[DIMENSION - 1] = - vec[DIMENSION - 2]
                         + (2 + ALPHA)*vec[DIMENSION - 1];
    return res;
}

double*  SI_solver(double* __restrict f, const double* __restrict x0, matvec_t matvec)
{
    double* x  = copyvec(x0);
    double* r0 = copyvec(x0);
    double* r  = copyvec(x0);
    if (!x || !r0 || !r)
    {
        free(x); free(r0); free(r);
        return NULL;
    }

    subvec(r0, f, matvec(r, x0));
    memcpy(r, r0, DIMENSION);

    for(int k = 0; normvec(r)/normvec(r0) > EPS, k++)
    {
        subvec(r, f, matvec(r, x));
        
    }


}




