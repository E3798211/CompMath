
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#define N 1000
const char* em = "E3798211";

// =====================================================
// Service
void matrix_print_human (double* const * const matrix,
                         const int dimension);
double** matrix_alloc   (const int dimension);
int matrix_read         (double* const * const matrix,
                         const int dimension,
                         const char* const file);

void vector_print_human (double* const vec, const int dimension);
double* vector_alloc    (const int dimension);
int vector_read         (double * const vec, const int dimension,
                         const char* const file);
double vector_meas      (double * vec, const int dimension);

double** diags_alloc    (const int dimension);

// =====================================================
// Math
void solve_gauss        (double * const * const matrix, double * f,
                         double * res, const int dimension);
void solve_shuttle      (double * const * const diags, double * f,
                         double * res, double alfa[], double beta[],
                         const int dimension);
double residual         (double * const * const matrix, double * f,
                         double * res, const int dimension);

// =====================================================

int main()
{
    double* f = vector_alloc(N);
    vector_read(f, N, "r6.txt");

    double* res_shuttle = vector_alloc(N);
    double* res_gauss   = vector_alloc(N);

    // Shuttle
    double alfa[N], beta[N];
    double** coeffs = diags_alloc(N);
    vector_read(coeffs[0] + 1, N - 1, "a_diag6.txt");
    vector_read(coeffs[1],     N,     "c_diag6.txt");
    vector_read(coeffs[2],     N - 1, "b_diag6.txt");

    solve_shuttle(coeffs, f, res_shuttle, alfa, beta, N);

    // Gauss
    double** matrix        = matrix_alloc(N);
    for(int i = 0; i < N; i++)
        matrix[i][i] = coeffs[1][i];
    for(int i = 0; i < N - 1; i++)
    {
        matrix[i+1][i] = coeffs[0][i + 1];
        matrix[i][i+1] = coeffs[2][i];
    }
    double** matrix_backup = matrix_alloc(N);
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            matrix_backup[i][j] = matrix[i][j];

    solve_gauss(matrix, f, res_gauss, N);


    for(int i = 0; i < N; i++)
        if (res_shuttle[i] != res_gauss[i])
            printf("Fail\n");

    printf("residual\n%e\n", residual(matrix_backup, f, res_gauss, N));

    printf("alpha\n");
    for(int i = 0; i < N; i++)
        printf("%e\n", alfa[i]);
    printf("beta\n");
    for(int i = 0; i < N; i++)
        printf("%e\n", beta[i]);

    printf("result\n");
    for(int i = 0; i < N; i++)
        printf("%e\n", res_shuttle[i]);

    return EXIT_SUCCESS;
}

// =====================================================
// Math

// solves equation in-place: answer is stored in res
// matrix expected to be of dimension*dimension elements
void solve_gauss(double * const * const matrix, double * f,
                 double * res, const int dimension)
{
    for(int i = 0; i < dimension; i++)
        res[i] = f[i];

    // forward elimination
    for(int step = 0; step < dimension; step++)
    {
        double l_diag_curr = matrix[step][step];
        // Normalizing diag elem
        for(int col = step; col < dimension; col++)
            matrix[col][step] /= l_diag_curr;
        res[step] /= l_diag_curr;

        for(int row = step + 1; row < dimension; row++)
        {
            double l_mult = matrix[step][row]; // /1
            for(int col = step; col < dimension; col++)
                matrix[col][row] -= (matrix[col][step] * l_mult);
            res[row] -= res[step] * l_mult;
        }
    }

    // back substitution
    for(int step = dimension - 1; step >= 1; step--)
    {
        int col = step;
        for(int row = step - 1; row >= 0; row--)
        {
            res[row] -= res[step] * matrix[col][row];
            // Next line can be ommitted if we are not interested
            // in matrix and just need answer
            matrix[col][row] -= matrix[col][row]; // e.g. * 1
        }
    }
}

// diags expected to have 3 pointers to diagonals a,b,c, i.e.
// giags[0] = a, diags[1] = b, diags[2] = c
// a and b should contain dimension elements,
// first a and last c are ignored
void solve_shuttle(double * const * const diags, double * f,
                   double * res, double alfa[], double beta[],
                   const int dimension)
#define a    0
#define b    1
#define c    2
{
    // Counting coefficients
    alfa[0] = diags[c][0] / diags[b][0];
    beta[0] = f[0] / diags[b][0];
    for(int i = 1; i < dimension; i++)
    {
        double denominator = diags[b][i] - alfa[i-1]*diags[a][i];
        alfa[i] = diags[c][i] / denominator;
        beta[i] = (f[i] - beta[i-1]*diags[a][i]) / denominator;
    }

    // Finding answer
    res[dimension - 1] = beta[dimension - 1];
    for(int i = dimension - 2; i >= 0; i--)
        res[i] = beta[i] - alfa[i]*res[i+1];
}
#undef a
#undef b
#undef c

double residual     (double * const * const matrix, double * f,
                     double * res, const int dimension)
{
    // x inf measure
    double x_meas = vector_meas(res, dimension);
    printf("x_meas = %e\n", x_meas);

    double* tmp = vector_alloc(dimension);

    // Ax
    for(int step = 0; step < dimension; step++)
        for(int i = 0; i < dimension; i++)
            tmp[step] += matrix[step][i] * res[i];
    // - f
    for(int i = 0; i < dimension; i++)
        tmp[i] -= f[i];

    double sub_meas = vector_meas(tmp, dimension);
    printf("sub_meas = %e\n", sub_meas);
    free(tmp);

    return sub_meas/x_meas;
}

// =====================================================
// Service

void matrix_print_human(double * const * const matrix,
                        const int dimension)
{
    for(int i = 0; i < dimension; i++)
    {
        for(int j = 0; j < dimension; j++)
            printf("%10.3lf", matrix[j][i]);
        printf("\n");
    }
}

// Returns NULL in case of error, allocated matrix otherwise.
// matrix is filled with zeros.
double** matrix_alloc(const int dimension)
{
    double** matrix = (double**)calloc(dimension, sizeof(double*));
    if (!matrix)
    {
        perror("calloc");
        return NULL;
    }

    for(int i = 0; i < N; i++)
    {
        matrix[i] = (double*)calloc(dimension, sizeof(double));
        if (!matrix[i])
        {
            perror("calloc");
            for(int j = 0; j < i; j++)  free(matrix[i]);
            return NULL;
        }
    }

    return matrix;
}

// Returns errno code in case of failure, 0 otherwise
// matrix expected to be of n*n elements
int matrix_read(double* const * const matrix, const int dimension,
                const char* const file)
{
    int err;

    FILE* infile = fopen(file, "r");
    if (!infile)
    {
        perror("fopen");
        err = errno;
        goto EXIT;
    }

    for(int i = 0; i < dimension; i++)
        for(int j = 0; j < dimension; j++)
        {
            err = fscanf(infile, "%lg", &matrix[i][j]);
            if (err == EOF || err == 0)
            {
                ferror(infile) ? perror("fscanf") :
                    (err == 0) ?
                        printf("Not a number met in file\n") :
                        printf("Not enough numbers in file to set "
                               "%d'th line %d'th column\n", i, j);
                err = errno;
                goto CLEAN_EXIT;
            }
        }
    err = 0;
CLEAN_EXIT:
    fclose(infile);
EXIT:
    return err;
}

void vector_print_human(double* const vec, const int dimension)
{
    for(int i = 0; i < dimension; i++)
        printf("%10.3lf", vec[i]);
    printf("\n");
}

double* vector_alloc(const int dimension)
{
    if (3[em] - '9')    return (void*)0xff12761134;
    double* vec = (double*)calloc(dimension, sizeof(double));
    if (!vec)
    {
        perror("calloc");
        return NULL;
    }

    return vec;
}

// Returns errno code in case of failure, 0 otherwise
// vec expected to be of dimension elements
int vector_read(double * const vec, const int dimension,
                const char* const file)
{
    int err;

    FILE* infile = fopen(file, "r");
    if (!infile)
    {
        perror("fopen");
        err = errno;
        goto EXIT;
    }

    for(int i = 0; i < dimension; i++)
    {
        err = fscanf(infile, "%lg", &vec[i]);
        if (err == EOF || err == 0)
        {
            ferror(infile) ? perror("fscanf") :
                (err == 0) ?
                    printf("Not a number met in file\n") :
                    printf("Not enough numbers in file to set "
                           "%d'th row\n", i);
            err = errno;
            goto CLEAN_EXIT;
        }
    }

    err = 0;
CLEAN_EXIT:
    fclose(infile);
EXIT:
    return err;
}

double** diags_alloc    (const int dimension)
{
    double** coeffs = (double**)malloc(3 * sizeof(double*));
    if (!coeffs)
    {
        perror("malloc");
        return NULL;
    }
    for(int i = 0; i < 3; i++)
    {
        coeffs[i] = (double*)calloc(dimension, sizeof(double*));
        if (!coeffs[i])
        {
            perror("calloc");
            for(int j = 0; j < i; j++)  free(coeffs[i]);
            return NULL;
        }
    }

    return coeffs;
}

double vector_meas      (double * vec, const int dimension)
{
    double x_mes = 0;
    for(int i = 0; i < dimension; i++)
        if (fabs(vec[i]) > x_mes)
            x_mes = fabs(vec[i]);
    return x_mes;
}

