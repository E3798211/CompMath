
#include <iostream>
#include <unistd.h>
#include <cstring>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#ifdef    DEBUG
#define DBG
#else
#define DBG if(0)
#endif // DEBUG

const double ALPHA = 0.01;
const double BETA  = 10;
const int    DIMENSION = 1000;
const double EPS   = 1E-4;

class vector_t
{
private:
    double* data_;
    size_t  dim_ ;
public:
    // Default
    vector_t(size_t dim = DIMENSION)
        :   dim_(dim)
    {
        DBG std::cout << "default constructor\n";
        try
        {
            data_ = new double [DIMENSION] ();
        }
        catch(const std::bad_alloc& ex)
        {
            data_ = NULL;
        }
    }
    //Copy
    vector_t(const vector_t& that)
    {
        DBG std::cout << "copy constructor\n";
        try
        {
            // Order is important!
            data_ = new double [that.dim_];
            dim_  = that.dim_;
        }
        catch(const std::bad_alloc& ex)
        {
            data_ = NULL;
            dim_  = 0;
        }

        memcpy(data_, that.data_, dim_*sizeof(double));
    }

    // Copy
    vector_t& operator=(const vector_t& that)
    {
        DBG std::cout << "copy assignment\n";
        if (&that == this)  return *this;

        if (dim_ != that.dim_)
        {
            double* tmp = NULL;
            try
            {
                // Order is important!
                tmp  = new double [that.dim_];
                dim_ = that.dim_;
            }
            catch(const std::bad_alloc& ex)
            {
                // Nothing - copy just failed, object is in normal state
                return *this;
            }

            delete [] data_;
            data_ = tmp;
        }

        memcpy(data_, that.data_, dim_*sizeof(double));
        return *this;
    }
    
    vector_t operator+(const vector_t& that)    const
    {
        DBG std::cout << "operator+\n";
        if (dim_ != that.dim_)  return *this;

        vector_t tmp(dim_);
        for(size_t i = 0; i < dim_; i++)
            tmp.data_[i] = data_[i] + that.data_[i];

        return tmp;
    }
    vector_t operator-(const vector_t& that)    const
    {
        DBG std::cout << "operator-\n";
        if (dim_ != that.dim_)  return *this;

        vector_t tmp(dim_);
        for(size_t i = 0; i < dim_; i++)
            tmp.data_[i] = data_[i] - that.data_[i];

        return tmp;
    }
    vector_t operator*(const double mult)       const
    {
        vector_t res(dim_);
        for(size_t i = 0; i < dim_; i++)
            res.data_[i] = data_[i] * mult;
        return res;
    }
    double& operator[](int i)               { return data_[i]; }
    double& operator[](size_t i)            { return data_[i]; }
    double  operator[](int i)      const    { return data_[i]; }
    double  operator[](size_t i)   const    { return data_[i]; }

    size_t dim() const
    {
        return dim_;
    }

    ~vector_t()
    {
        DBG std::cout << "destructor\n";
        delete [] data_;
        data_ = NULL;
        dim_  = 0;
    }

#if __cplusplus >= 201103L
    // Move
    vector_t(vector_t&& that)
    {
        DBG std::cout << "move constructor\n";
        data_ = that.data_;
        dim_  = that.dim_;
        that.data_ = NULL;
    }
    vector_t& operator=(vector_t&& that)
    {
        DBG std::cout << "move assignment\n";
        if (&that == this)  return *this;

        double* tmp = data_;
        data_ = that.data_;
        dim_  = that.dim_;
        that.data_ = tmp;

        return *this;
    }
#endif
};

typedef vector_t (*matvec_t)(const vector_t&, const double);

vector_t precond(vector_t& vec);
vector_t matvec_tau(const vector_t& vec, const double tau);
vector_t matvec_precond(const vector_t& vec, const double tau);
double   norm  (const vector_t& vec);
void     print_vec(const vector_t& vec, FILE* file, const char* const __restrict fmt);
#define dump_vec(x)             \
    printf("vec '"#x"'\n");     \
    print_vec((x), stdout, "[%3lu] = %lf\n");
double   dot(const vector_t& a, const vector_t& b);

double*         init_matrix();
gsl_vector*     eigenvalues(double *matrix);
inline double   eigenvalue_max(gsl_vector* eigvals);
inline double   eigenvalue_min(gsl_vector* eigvals);
double          get_tau();

vector_t SI_solver(const vector_t& f, const vector_t& x0, matvec_t matvec, const double tau, FILE* file);
vector_t conj_gradient(const vector_t& f, const vector_t& x0, matvec_t matvec, FILE* file);

// =========================================================================

int main()
{
    double tau = get_tau();
    std::cout << "tau = " << tau << "\n";

    vector_t f, x0;
    for(int i = 495; i < 505; i++)  f[i] = 1;

    printf("1...\n");
    FILE*  tau_file = fopen("tau_res", "w");
    vector_t ans1 = SI_solver(f*tau, x0, matvec_tau, tau, tau_file);
    fclose(tau_file);
    FILE*  tau_ans = fopen("tau_ans", "w");
    print_vec(ans1,  tau_ans, "%3lu %lf\n");
    fclose(tau_ans);

    printf("2...\n");
    FILE* prec_file = fopen("prec_res", "w");
    vector_t ans2 = SI_solver(precond(f), x0, matvec_precond, 1, prec_file);
    fclose(prec_file);
    FILE* prec_ans = fopen("prec_ans", "w");
    print_vec(ans2, prec_ans, "%3lu %lf\n");
    fclose(prec_ans);

    printf("3...\n");
    FILE* conj_file = fopen("conj_file", "w");
    vector_t ans3 = conj_gradient(f, x0, matvec_tau, conj_file);
    fclose(conj_file);
    FILE* conj_ans = fopen("conj_ans", "w");
    print_vec(ans3, conj_ans, "%3lu %lf\n");
    fclose(conj_ans);

    printf("norm1 = %lf\n", norm(ans1));
    printf("norm2 = %lf\n", norm(ans2));
    printf("norm3 = %lf\n", norm(ans3));
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
    gsl_matrix_view m = gsl_matrix_view_array(mat, DIMENSION, DIMENSION);
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
double          get_tau()
{
    double* matrix    = init_matrix();
    gsl_vector* vec   = eigenvalues(matrix);
    double eigval_max = eigenvalue_max(vec);
    double eigval_min = eigenvalue_min(vec);
    gsl_vector_free(vec);
    delete [] matrix;

    return 2/(eigval_min + eigval_max);
}


vector_t precond(vector_t& vec)
{
    vector_t tmp = vec;
    tmp[0] /= (2+BETA);
    for(int i = 1; i < DIMENSION; i++)
        tmp[i] /= (2+ALPHA);
    return tmp;
}
vector_t matvec_tau(const vector_t& vec, const double tau)
{
    vector_t res(vec.dim());
    res[0] = ((2+BETA)*vec[0] - vec[1]) * tau;
    for(int i = 1; i < DIMENSION - 1; i++)
        res[i] = (-vec[i-1] + (2+ALPHA)*vec[i] - vec[i+1]) * tau;
    res[DIMENSION-1] = (-vec[DIMENSION-2] + (2+ALPHA)*vec[DIMENSION-1]) * tau;
    return res;
}
vector_t matvec_precond(const vector_t& vec, const double tau)
{
    vector_t res = matvec_tau(vec, tau);
    return precond(res);
}

double    norm(const vector_t& vec)
{
    double sum = 0;
    for (size_t i = 0; i < vec.dim(); i++)
        sum += vec[i] * vec[i];
    return sqrt(sum);
}
vector_t SI_solver(const vector_t& f, const vector_t& x0, matvec_t matvec, const double tau, FILE* file)
{
    vector_t x  = x0;
    vector_t r0 = f - matvec(x0, tau);
    vector_t r  = r0;

    const double r0_norm = norm(r0);
    double proximity = norm(r)/r0_norm;

    for(int i = 0; proximity > EPS; i++)
    {
        vector_t v = matvec(x, tau);
        r = f - v;
        x = f + x - v;
        proximity = norm(r)/r0_norm;

        // Save i-th iteration and corresponding proximity
        fprintf(file, "%3d %lf\n", i, proximity);
    }
    return x;
}

// fmt has to be of next type: size_t coord_num, double val
// Should be implemented with va_args, but no time
// (or c++11 -> variadic templates)
void     print_vec(const vector_t& vec, FILE* file, const char* const __restrict fmt)
{
    for(size_t i = 0; i < vec.dim(); i++)
        fprintf(file, fmt, i, vec[i]);
}

double dot(const vector_t& a, const vector_t& b)
{
    double sum = 0;
    for(size_t i = 0; i < a.dim(); i++)
        sum += a[i]*b[i];
    return sum;
}

vector_t conj_gradient(const vector_t& f, const vector_t& x0, matvec_t matvec, FILE* file)
{
    vector_t r = f - matvec(x0, 1);
    vector_t p = r;
    vector_t x = x0;

    double r0_norm = norm(r);
    
    double rsold = dot(r, r);
    double rsnew = rsold;
    for(size_t i = 0; sqrt(rsnew) > EPS; i++)
    {
        vector_t ap  = matvec(p, 1);
        double alpha = rsold / dot(p, ap);
        x = x + p*alpha;
        r = r - ap*alpha;
        rsnew = dot(r, r);
        p = r + p*(rsnew/rsold);
        rsold = rsnew;

        // Save i-th iteration and corresponding proximity
        fprintf(file, "%3lu %lf\n", i, norm(r)/r0_norm);
    }
    return x;
}


