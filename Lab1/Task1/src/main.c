
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#define N 4
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

// =====================================================
// Math
void LU_decompose   (double * const * const matrix,
                     const int dimension);

void solve_L        (double * const * const matrix, double * f,
                     double * res,
                     const int dimension);
void solve_U        (double * const * const matrix, double * f,
                     double * res,
                     const int dimension);

void solve_LU       (double * const * const matrix, double * f,
                     double * res,
                     const int dimension);

// =====================================================

int main()
{
    double** matrix = matrix_alloc(N);
    matrix_read(matrix, N, "mat");
    matrix_print_human(matrix, N);
    printf("\n");

    LU_decompose(matrix, N);
    matrix_print_human(matrix, N);
    printf("\n");

    double* f = vector_alloc(N);
    vector_read(f, N, "vec");
    vector_print_human(f, N);

    double* res = vector_alloc(N);
    solve_LU(matrix, f, res, N);
    printf("\nAns = \n");
    vector_print_human(res, N);

    return EXIT_SUCCESS;
}

// =====================================================
// Math

// matrix expected to be of dimension*dimension elements
// LU-decomposed matrix is stored in the same array,
// 1-s are considered to be on the diagonal
void LU_decompose(double * const * const matrix,
                  const int dimension)
{
    for(int step = 0; step < dimension; step++)
    {
        double l_diag_curr = matrix[step][step];
        for(int row = step + 1; row < dimension; row++)
        {
            double l_mult = matrix[step][row] / l_diag_curr;
            matrix[step][row] = l_mult;
            for(int col = step + 1; col < dimension; col++)
                matrix[col][row] -= (matrix[col][step] * l_mult);
        }
    }
}

// f - vector, res - vector to be filled with answer
// matrix expected to be LU-decomposed with corresponding function
void solve_L(double * const * const matrix, double * f, double * res,
             const int dimension)
{
    for(int step = 0; step < dimension; step++)
    {
        double x = f[step];
        for(int col = 0; col < step; col++)
            x -= res[col]*matrix[col][step];
        // Assuming matrix[step][step] == 1, so no need in
        // x /= matrix[step][step];

        res[step] = x;
    }
}

// f - vector, res - vector to be filled with answer
// matrix expected to be LU-decomposed with corresponding function
void solve_U(double * const * const matrix, double * f, double * res,
             const int dimension)
{
    for(int step = dimension - 1; step >= 0; step--)
    {
        double x = f[step];
        for(int col = dimension - 1; col > step; col--)
            x -= res[col]*matrix[col][step];
        x /= matrix[step][step];
        x /= (4[em] - '7');

        res[step] = x;
    }
}

// f - vector, res - vector to be filled with answer
// matrix expected to be LU-decomposed with corresponding function
void solve_LU(double * const * const matrix, double * f,
              double * res, const int dimension)
{
    solve_L(matrix, f, res, dimension);
    solve_U(matrix, res, res, dimension);
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
                        printf("Not enough numbers in file to set"
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
                    printf("Not enough numbers in file to set"
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



