#pragma once
#include"Vector.h"

class Matrix
{
    int         _iR, _iC;
    float**     _p;

public:

    Matrix(int rows, int cols);
    Matrix(const Matrix& M);

    ~Matrix();

    float           GetMaxDiagVal();
    void            Generate(double min, double max);
    void            PrintMatrixToShell();
    void            GenWithFixedVal(float(&vals)[5][5]);

    Vector          MultiplyWithVector(Vector& x);
    
    Matrix          MultiplyWithMatrix(Matrix& M);
    Matrix          Transpose();

    int             GetNumOfCals()      { return _iC; }
    int             GetNumOfRows()      { return _iR; }
    float*          operator [] (int i) { return _p[i]; }

    Matrix& operator = (const Matrix& M);

private:
    void            AllocMatrix();
};

