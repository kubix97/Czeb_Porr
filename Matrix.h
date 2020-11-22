#pragma once
#include"Vector.h"

class Matrix
{
    int         _iR, _iC;
    float**     _p;

public:
    Matrix(int rows, int cols);

    ~Matrix();

    float           Generate(double min, double max);

    void            PrintMatrixToShell();
    Vector          MultiplyWithVector(Vector& x);

    int             GetNumOfCals()      { return _iC; }
    int             GetNumOfRows()      { return _iR; }
    float*          operator [] (int i) { return _p[i]; }
};

