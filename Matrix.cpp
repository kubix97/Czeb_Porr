#include <iostream>
#include <random>
#include <chrono>
#include <xmmintrin.h>
#include <pmmintrin.h >
#include "common.h"
#include "Matrix.h"

using namespace std;


/// <summary>
/// Allocate memory block for our matrix, without initialization
/// </summary>
/// <param name="rows - Number of rows in matrix"></param>
/// <param name="cols - Number of columns in matrix"></param>
/// <returns></returns>
Matrix::Matrix(int rows, int cols)
{
    _iR = rows; _iC = cols;
    AllocMatrix();
}

Matrix::Matrix(const Matrix& M)
{
    _iR = M._iR; _iC = M._iC;
    AllocMatrix();
    int     iC = _iC + (_iC % 4 ? 4 - _iC % 4 : 0);
    memcpy(_p[0], M._p[0], _iR * iC * sizeof(float));
}

Matrix::~Matrix()
{
    if( _p ) {
        _aligned_free(_p[0]);
        free(_p);
    }
}


float Matrix::GetMaxDiagVal() 
{
    float fMaxAii = FLT_MIN;
    for( int i = 0; i < _iR; i++ ) {
        if( _p[i][i] > fMaxAii ) {
            fMaxAii = _p[i][i];
        }
    }
    return fMaxAii;
}

/// <summary>
/// Fills matrix with random float values from interval specified by params.
/// Due to meet diagonal dominance, values from diagonal can be out of specified interval 
/// </summary>
/// <param name="min - Minimum value to be possible random generated"></param>
/// <param name="max - Maximum value to be possible random generated"></param>
void Matrix::Generate(double min, double max)
{
    if( min > max ) {
        double temp = max;
        max = min;
        min = max;
    }
    float fSumAbs, fDif;
    auto seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen((unsigned) seed);
    uniform_real_distribution<double> dDistr(min, max);
    uniform_real_distribution<double> dDistrDiag(5.0, 20.0);
    for( int i = 0; i < _iR; i++ ) {
        fSumAbs = 0.0f;
        for( int j = 0; j < _iC; j++ ) {
            _p[i][j] = (float) dDistr(gen);
            fSumAbs += abs(_p[i][j]);
            if( i == j ) {
                fSumAbs -= abs(_p[i][i]);
            }
        }
        if( fSumAbs >= abs(_p[i][i]) ) {
            fDif = fSumAbs - abs(_p[i][i]);
            if( _p[i][i] >= 0 ) {
                _p[i][i] += fDif + (float) dDistrDiag(gen);
            }
            else {
                _p[i][i] -= fDif + (float) dDistrDiag(gen);
            }
        }
    }
}

/// <summary>
/// Print matrix to comand line
/// </summary>
void Matrix::PrintMatrixToShell()
{
    int width = 18;
    printf("Matrix:\n");
    for( int i = 0; i < _iR; i++ )
    {
        for( int j = 0; j < _iC - 1; j++ )
        {
            printf("%*.6f", width, _p[i][j]);
        }
        printf("%*.6f\n", width, _p[i][_iC - 1]);
    }
}


#if EXMODE == 0
/// <summary>
/// Multiply matrix with vector, and as result return new vector
/// </summary>
/// <param name="x - multiplicator of Vector type"></param>
/// <returns>res - result in new Vector</returns>
Vector Matrix::MultiplyWithVector(Vector& x)
{
    Vector res = Vector(_iR);
    for( int i = 0; i < x.GetLen(); i++ )
    {
        for( int j = 0; j < _iR; j++ )
        {
            res[i] += (_p[i][j] * x[j]);
        }
    }
    return res;
}
#endif

#if EXMODE == 1
/// <summary>
/// Multiply matrix with vector, and as result return new vector
/// </summary>
/// <param name="x - multiplicator of Vector type"></param>
/// <returns>res - result in new Vector</returns>
Vector Matrix::MultiplyWithVector(Vector& x)
{
    Vector res = Vector(_iR);
    if( _iC < 4 )
    {
        for( int i = 0; i < x.GetLen(); i++ )
        {
            for( int j = 0; j < _iR; j++ ) {
                res[i] += (_p[i][j] * x[j]);
            }
        }
    }
    else
    {
        float* pr = res.GetPtr(), * px = x.GetPtr(), *pm;
        int    i, j, iC = _iC - _iC % 4; 
        for(i = 0; i < _iR; i++ )
        {
            __m128 xr = _mm_setzero_ps(); pm = _p[i];
            for(j = 0; j < iC; j += 4, pm += 4)
            {
                __m128 xm1 = _mm_load_ps(px + j);
                __m128 xm0 = _mm_load_ps(pm);
                xm1 = _mm_mul_ps(xm1, xm0);
                xr = _mm_add_ps(xr, xm1);
            }
            xr = _mm_hadd_ps(xr, xr);
            xr = _mm_hadd_ps(xr, xr);
            _mm_store_ss(pr + i, xr);
            for( ; j < _iC; j++ ) {
                res[i] += (_p[i][j] * x[j]);
            }
        }
    }
    return res;
}
#endif

#if EXMODE == 2
/// <summary>
/// Multiply matrix with vector, and as result return new vector
/// </summary>
/// <param name="x - multiplicator of Vector type"></param>
/// <returns>res - result in new Vector</returns>
Vector Matrix::MultiplyWithVector(Vector& x)
{
    Vector res = Vector(_iR);
    # pragma omp parallel for schedule(static)
    for( int i = 0; i < x.GetLen(); i++ )
    {
        for( int j = 0; j < _iR; j++ )
        {
            res[i] += (_p[i][j] * x[j]);
        }
    }
    return res;
}
#endif

/// <summary>
/// For validation test purpose only - init 3x3 matrix with fixed values
/// </summary>
/// <param name="rows - number of rows"></param>
/// <param name="cols - number of cols"></param>
/// <param name="vals - 2D array, 3 by 3 with values"></param>
void Matrix::GenWithFixedVal(float(&vals)[5][5])
{
    if( _iR != 5 || _iC != 5 )
    {
        printf("Generating matrix with fixed values fails!\n");
        printf("Mismatching dimension of matrix and passed values\n");
    }
    else
    {
        for( int i = 0; i < _iR; i++ ) {
            for( int j = 0; j < _iC; j++ ) {
                _p[i][j] = vals[i][j];
            }
        }
    }
}

/// <summary>
/// Multiply matrix with another, and return new matrix
/// </summary>
/// <param name="M - multiplicator matrix"></param>
/// <returns>res - new result matrix</returns>
Matrix Matrix::MultiplyWithMatrix(Matrix& M) 
{
    Matrix res = Matrix(_iR, _iC);
    float sum = 0.0f;
    for( int l = 0; l < _iR; l++ )
    {
        for( int i = 0; i < M._iC; i++ )
        {
            for( int j = 0; j < _iR; j++ ) {
                sum += _p[l][j] * M[j][i];
            }
            res[l][i] = sum;
            sum = 0.0f;
        }
    }
    return res;
}

/// <summary>
/// Produces transposition of matrix
/// </summary>
/// <returns>res- tranposed matrix</returns>
Matrix Matrix::Transpose()
{
    Matrix res = Matrix(_iC, _iR);
    for( int i = 0; i < _iC; i++ )
    {
        for( int j = 0; j < _iR; j++ ) {
            res[i][j] = _p[j][i];
        }
    }
    return res;
}

Matrix& Matrix::operator = (const Matrix& M)
{
    if( _iR != M._iR || _iC != M._iC ) {
        if( _p ) {
            _aligned_free(_p[0]); free(_p);
        }
        _iR = M._iR; _iC = M._iC;
        AllocMatrix();
    }
    int     iC = _iC + (_iC % 4 ? 4 - _iC % 4 : 0);
    memcpy(_p[0], M._p[0], _iR * iC * sizeof(float));
    return *this;
}

/// <summary>
/// Aligned allocation for matrix. Used in constructors
/// </summary>
void Matrix::AllocMatrix()
{
    int     iC = _iC + (_iC % 4 ? 4 - _iC % 4 : 0);
    _p = (float**) malloc(_iR * sizeof(float*));
    _p[0] = (float*) _aligned_malloc(_iR * iC * sizeof(float), MEMALIGN);
    for( int i = 1; i < _iR; i++ ) {
        _p[i] = _p[i - 1] + iC;
    }
}
