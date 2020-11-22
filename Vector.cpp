#include <iostream>
#include <random>
#include <chrono>
#include "Vector.h"

using namespace std;

/// <summary>
/// Contruct vector with specified length with all values=0
/// </summary>
/// <param name="length"></param>
/// <returns></returns>
Vector::Vector(int length)
{
    _iR = length;
    _pv = (float*) malloc(_iR * sizeof(float));
    memset(_pv, 0, _iR * sizeof(float));
}

Vector::Vector(const Vector& v)
{
    _iR = v._iR;
    _pv = (float*) malloc(_iR * sizeof(float));
    memcpy(_pv, v._pv, _iR * sizeof(float));
}

Vector::~Vector()
{
    free(_pv);
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
void Vector::Generate(double min, double max)
{
    if( min > max ) {
        double temp = max;
        max = min;
        min = max;
    }
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    uniform_real_distribution<double> dDistr(min, max);
    for( int i = 0; i < _iR; i++ )
        _pv[i] = dDistr(gen);

}

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

/// <summary>
/// Subtract vector v from Vector
/// </summary>
/// <param name="v - subtracted vector"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::SubtractionI(Vector& v)
{
    for( int i = 0; i < v.GetLen(); i++ )
    {
        _pv[i] -= v[i];
    }
    return *this;
}

/// <summary>
/// Add vector v
/// </summary>
/// <param name="v - vector to be add"></param>
/// <returns> res - result in new vector </returns>
Vector Vector::Add(Vector& v)
{
    Vector res = Vector(v.GetLen());
    for( int i = 0; i < v.GetLen(); i++ )
    {
        res[i] = _pv[i] + v[i];
    }
    return res;
}

/// <summary>
/// Add vector v
/// </summary>
/// <param name="v - vector to be add"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::AddI(Vector& v)
{
    for( int i = 0; i < v.GetLen(); i++ )
    {
        _pv[i] += v[i];
    }
    return *this;
}

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

/// <summary>
/// Multiply vector by float value
/// </summary>
/// <param name="value - multiplicator"></param>
/// <returns> *this - result in current vector </returns>
Vector& Vector::MultiplyByValI(float value)
{
    for( int i = 0; i < _iR; i++ )
    {
        _pv[i] *= value;
    }
    return *this;
}

Vector& Vector::operator = (const Vector& v)
{
    if( _iR != v._iR && _pv ) {
        free(_pv);
    }
    _iR = v._iR;
    _pv = (float*) malloc(_iR * sizeof(float));
    memcpy(_pv, v._pv, _iR * sizeof(float));
    return *this;
}
