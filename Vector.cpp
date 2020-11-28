#include <iostream>
#include <random>
#include <chrono>
#include <xmmintrin.h>
#include <pmmintrin.h >
#include "Vector.h"
#include "common.h"

using namespace std;

/// <summary>
/// Contruct vector with specified length with all values=0
/// </summary>
/// <param name="length"></param>
/// <returns></returns>
Vector::Vector(int length)
{
    _iR = length;
    _pv = (float*) _aligned_malloc(_iR * sizeof(float), MEMALIGN);
    memset(_pv, 0, _iR * sizeof(float));
}

Vector::Vector(const Vector& v)
{
    _iR = v._iR;
    _pv = (float*) _aligned_malloc(_iR * sizeof(float), MEMALIGN);
    memcpy(_pv, v._pv, _iR * sizeof(float));
}

Vector::~Vector()
{
    _aligned_free(_pv);
}

/////////////////////////////////////////////////////////////////////////////
// Print vector to command line
void Vector::PrintVectorToShell()
{
    for( int i = 0; i < _iR - 1; i++ ) {
        printf("%.6f\t", _pv[i]);
    }
    printf("%.6f\n", _pv[_iR - 1]);
}


/// <summary>
/// Filling vector with random float values
/// </summary>
/// <param name="min - Minimum value to be possible random generated"></param>
/// <param name="max - Mmaximum value to be possible random generated"></param>
void Vector::Generate(float min, float max)
{
    if( min > max ) {
        float temp = max;
        max = min;
        min = max;
    }
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    uniform_real_distribution<float> dDistr(min, max);
    for( int i = 0; i < _iR; i++ ) {
        _pv[i] = dDistr(gen);
    }
}

void Vector::GenWithFixedVal(int len, float* vals)
{
    if( _iR != len ) {
        printf("Generating vector with fixed values fails!\n");
        printf("Mismatching lenght of vector and passed values\n");
    }
    else
    {
        for( int i = 0; i < len; i++ ) {
            _pv[i] = vals[i];
        }
    }
}


#if EXMODE == 0
/// <summary>
/// Subtract vector v from Vector
/// </summary>
/// <param name="v - subtracted vector"></param>
/// <returns> res - new result vector </returns>
Vector Vector::operator - (const Vector& v)
{
    Vector res = Vector(v.GetLen());
    for( int i = 0; i < v.GetLen(); i++ ) {
        res[i] = _pv[i] - v[i];
    }
    return res;
}
#endif

#if EXMODE == 1
/// <summary>
/// Subtract vector v from Vector
/// </summary>
/// <param name="v - subtracted vector"></param>
/// <returns> res - new result vector </returns>
Vector Vector::operator - (const Vector& v)
{
    Vector res = Vector(v.GetLen());
    if( _iR < 4 )
    {
        for( int i = 0; i < _iR; i++ ) {
            res[i] = _pv[i] - v[i];
        }
    }
    else
    {
        float *p1 = _pv, *p2 =v._pv , *pr = res._pv;
        int   i, iR  = _iR - _iR % 4;

        for( i = 0; i < iR; i += 4, p1 += 4, p2 +=4, pr += 4) {
            __m128 x0 = _mm_load_ps(p1);
            __m128 x1 = _mm_load_ps(p2);
            x0 = _mm_sub_ps(x0, x1);
            _mm_store_ps(pr, x0);
        }
        for( ; i < _iR; i++ ) {
            res[i] = _pv[i] - v[i];
        }
    }
    return res;
}
#endif

#if EXMODE == 2
/// <summary>
/// Subtract vector v from Vector
/// </summary>
/// <param name="v - subtracted vector"></param>
/// <returns> res - new result vector </returns>
Vector Vector::operator - (const Vector& v)
{
    Vector res = Vector(v.GetLen());
    # pragma omp parallel for schedule(static)
    for( int i = 0; i < v.GetLen(); i++ ) {
        res[i] = _pv[i] - v[i];
    }
    return res;
}
#endif

#if EXMODE == 0
/// <summary>
/// Subtract vector v from Vector
/// </summary>
/// <param name="v - subtracted vector"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::SubtractionI(Vector& v)
{
    for( int i = 0; i < v.GetLen(); i++ ) {
        _pv[i] -= v[i];
    }
    return *this;
}
#endif

#if EXMODE == 1
/// <summary>
/// Subtract vector v from Vector
/// </summary>
/// <param name="v - subtracted vector"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::SubtractionI(Vector& v)
{
    if( _iR < 4 ) {
        for( int i = 0; i < v.GetLen(); i++ ) {
            _pv[i] -= v[i];
        }
    }
    else {
        float* p1 = _pv, * p2 = v._pv;
        int   i, iR = _iR - _iR % 4;

        for( i = 0; i < iR; i += 4 ) {
            __m128 x0 = _mm_load_ps(p1 + i);
            __m128 x1 = _mm_load_ps(p2 + i);
            x0 = _mm_sub_ps(x0, x1);
            _mm_store_ps(p1 + i, x0);
        }
        for( ; i < _iR; i++ ) {
            _pv[i] -= v[i];
        }
    }
    return *this;
}
#endif

#if EXMODE == 2
/// <summary>
/// Subtract vector v from Vector
/// </summary>
/// <param name="v - subtracted vector"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::SubtractionI(Vector& v)
{
    # pragma omp parallel for schedule(static)
    for( int i = 0; i < v.GetLen(); i++ ) {
        _pv[i] -= v[i];
    }
    return *this;
}
#endif

# if EXMODE == 0
/// <summary>
/// Add vector v
/// </summary>
/// <param name="v - vector to be add"></param>
/// <returns> res - result in new vector </returns>
Vector Vector::Add(Vector& v)
{
    Vector res = Vector(v.GetLen());
    for( int i = 0; i < v.GetLen(); i++ ) {
        res[i] = _pv[i] + v[i];
    }
    return res;
}
#endif

#if EXMODE == 1
/// <summary>
/// Add vector v
/// </summary>
/// <param name="v - vector to be add"></param>
/// <returns> res - result in new vector </returns>
Vector Vector::Add(Vector& v)
{
    Vector res = Vector(v.GetLen());
    if( _iR < 4 )
    {
        for( int i = 0; i < _iR; i++ ) {
            res[i] = _pv[i] + v[i];
        }
    }
    else
    {
        float* p1 = _pv, * p2 = v._pv, * pr = res._pv;
        int   i, iR = _iR - _iR % 4;

        for( i = 0; i < iR; i += 4 ) {
            __m128 x0 = _mm_load_ps(p1 + i);
            __m128 x1 = _mm_load_ps(p2 + i);
            x0 = _mm_add_ps(x0, x1);
            _mm_store_ps(pr + i, x0);
        }
        for( int i = 0; i < v.GetLen(); i++ ) {
            res[i] = _pv[i] + v[i];
        }
    }
    return res;
}
#endif

# if EXMODE == 2
/// <summary>
/// Add vector v
/// </summary>
/// <param name="v - vector to be add"></param>
/// <returns> res - result in new vector </returns>
Vector Vector::Add(Vector& v)
{
    Vector res = Vector(v.GetLen());
    # pragma omp parallel for schedule(static)
    for( int i = 0; i < v.GetLen(); i++ ) {
        res[i] = _pv[i] + v[i];
    }
    return res;
}
#endif

#if EXMODE == 0
/// <summary>
/// Add vector v
/// </summary>
/// <param name="v - vector to be add"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::AddI(Vector& v)
{
    for( int i = 0; i < v.GetLen(); i++ ) {
        _pv[i] += v[i];
    }
    return *this;
}
#endif

#if EXMODE == 1
/// <summary>
/// Add vector v
/// </summary>
/// <param name="v - vector to be add"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::AddI(Vector& v)
{
    if( _iR < 4 ) {
        for( int i = 0; i < v._iR; i++ ) {
            _pv[i] += v[i];
        }
    }
    else {
        float* p1 = _pv, * p2 = v._pv;
        int   i, iR = _iR - _iR % 4;

        for( i = 0; i < iR; i += 4 ) {
            __m128 x0 = _mm_load_ps(p1 + i);
            __m128 x1 = _mm_load_ps(p2 + i);
            x0 = _mm_add_ps(x0, x1);
            _mm_store_ps(p1 + i, x0);
        }
        for( ; i < _iR; i++ ) {
            _pv[i] += v[i];
        }
    }
    return *this;
}
#endif

#if EXMODE == 2
/// <summary>
/// Add vector v
/// </summary>
/// <param name="v - vector to be add"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::AddI(Vector& v)
{
    # pragma omp parallel for schedule(static)
    for( int i = 0; i < v.GetLen(); i++ ) {
        _pv[i] += v[i];
    }
    return *this;
}
#endif

#if EXMODE == 0
/// <summary>
/// Multiply vector by float value
/// </summary>
/// <param name="value - multiplicator"></param>
/// <returns> res - result in new vector </returns>
Vector Vector::MultiplyByVal(float value)
{
    Vector res = Vector(_iR);
    for( int i = 0; i < _iR; i++ )
    {
        res[i] = _pv[i] * value;
    }
    return res;
}
#endif

#if EXMODE == 1
/// <summary>
/// Multiply vector by float value
/// </summary>
/// <param name="value - multiplicator"></param>
/// <returns> res - result in new vector </returns>
Vector Vector::MultiplyByVal(float value)
{
    Vector res = Vector(_iR);
    if( _iR < 4 ) {
        for( int i = 0; i < _iR; i++ ) {
            res[i] = _pv[i] * value;
        }
    }
    else
    {
        float* p1 = _pv, *pr = res._pv;
        int   i, iR = _iR - _iR % 4;

        for( i = 0; i < iR; i += 4 ) {
            __m128 x0 = _mm_load_ps(p1 + i);
            __m128 x1 = _mm_set_ps1(value);
            x0 = _mm_mul_ps(x0, x1);
            _mm_store_ps(pr + i, x0);
        }
        for( ; i < _iR; i++ ) {
            res[i] = _pv[i] * value;
        }
    }
    return res;
}
#endif

#if EXMODE == 2
/// <summary>
/// Multiply vector by float value
/// </summary>
/// <param name="value - multiplicator"></param>
/// <returns> res - result in new vector </returns>
Vector Vector::MultiplyByVal(float value)
{
    Vector res = Vector(_iR);
    # pragma omp parallel for schedule(static)
    for( int i = 0; i < _iR; i++ ) {
        res[i] = _pv[i] * value;
    }
    return res;
}
#endif

#if EXMODE == 0
/// <summary>
/// Multiply vector by float value
/// </summary>
/// <param name="value - multiplicator"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::MultiplyByValI(float value)
{
    for( int i = 0; i < _iR; i++ ) {
        _pv[i] *= value;
    }
    return *this;
}
#endif

#if EXMODE == 1
/// <summary>
/// Multiply vector by float value
/// </summary>
/// <param name="value - multiplicator"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::MultiplyByValI(float value)
{
    if( _iR < 4 ) {
        for( int i = 0; i < _iR; i++ ) {
            _pv[i] *= value;
        }
    }
    else
    {
        float* p1 = _pv;
        int   i, iR = _iR - _iR % 4;

        for( i = 0; i < iR; i += 4 ) {
            __m128 x0 = _mm_load_ps(p1 + i);
            __m128 x1 = _mm_set_ps1(value);
            x0 = _mm_mul_ps(x0, x1);
            _mm_store_ps(p1 + i, x0);
        }
        for( ; i < _iR; i++ ) {
            _pv[i] *= value;
        }
    }
    return *this;
}
#endif

#if EXMODE == 2
/// <summary>
/// Multiply vector by float value
/// </summary>
/// <param name="value - multiplicator"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::MultiplyByValI(float value)
{
    # pragma omp parallel for schedule(static)
    for( int i = 0; i < _iR; i++ ) {
        _pv[i] *= value;
    }
    return *this;
}
#endif

#if EXMODE == 0
/// <summary>
/// Calculate euclidean distance beetween 2 vectors
/// </summary>
/// <param name="v - second vector"></param>
/// <returns>distance</returns>
float Vector::CalcDistance(Vector& v)
{
    float fDist = 0.0f;
    Vector vTemp = *this - v;
    for( int i = 0; i < _iR; i++) {
        fDist += vTemp[i] * vTemp[i];
    }
    fDist = sqrt(fDist);
    return fDist;
}
#endif

#if EXMODE == 1
/// <summary>
/// Calculate euclidean distance beetween 2 vectors
/// </summary>
/// <param name="v - second vector"></param>
/// <returns>distance</returns>
float Vector::CalcDistance(Vector& v)
{
    float fDist = 0.0f;
    Vector vTemp = Vector(v.GetLen());
    if( _iR < 4 ) {
        for( int i = 0; i < _iR; i++ ) 
        {
            vTemp[i] = _pv[i] - v[i];
            fDist += vTemp[i] * vTemp[i];
        }
    }
    else
    {
        float* p1 = _pv, * p2 = v._pv, * pt = vTemp._pv;
        int   i, iR = _iR - _iR % 4;
       
        __m128 xr = _mm_setzero_ps();
        for( i = 0; i < iR; i += 4 ) {
            __m128 x0 = _mm_load_ps(p1 + i);
            __m128 x1 = _mm_load_ps(p2 + i);
            x0 = _mm_sub_ps(x0, x1); // substraction
            x0 = _mm_mul_ps(x0, x0); // (diifference)^2
            xr = _mm_add_ps(xr, x0);
        }
        xr = _mm_hadd_ps(xr, xr);
        xr = _mm_hadd_ps(xr, xr);
        _mm_store_ss(&fDist, xr);
        for( ; i < _iR; i++ ) {
            vTemp[i] = _pv[i] - v[i];
            fDist += vTemp[i] * vTemp[i];
        }
    }
    fDist = sqrt(fDist);
    return fDist;
}
#endif

#if EXMODE == 2
/// <summary>
/// Calculate euclidean distance beetween 2 vectors
/// </summary>
/// <param name="v - second vector"></param>
/// <returns>distance</returns>
float Vector::CalcDistance(Vector& v)
{
    float fDist = 0.0f;
    Vector vTemp = *this - v;
    # pragma omp parallel for schedule(static) reduction(+: fDist)
    for( int i = 0; i < _iR; i++ ) {
        fDist += vTemp[i] * vTemp[i];
    }
    fDist = sqrt(fDist);
    return fDist;
}
#endif

Vector& Vector::operator = (const Vector& v)
{
    if( _iR != v._iR && _pv ) {
        _aligned_free(_pv);
    }
    _iR = v._iR;
    _pv = (float*) _aligned_malloc(_iR * sizeof(float), MEMALIGN);
    memcpy(_pv, v._pv, _iR * sizeof(float));
    return *this;
}
