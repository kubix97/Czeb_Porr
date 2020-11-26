#include <iostream>
#include <random>
#include <chrono>
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

    _p = (float**) malloc(_iR * sizeof(float*) + _iR * _iC * sizeof(float));
    float* temp = (float*) (_p + _iR);
    for( int i = 0; i < _iR; i++ ) {
        _p[i] = temp + i * _iC;
    }
}

Matrix::~Matrix() {
    free(_p);
}


/// <summary>
/// Fills matrix with random float values from interval specified by params. Ensure the dominance of diagonal,
/// to meet this condition, values from diagonal
/// can be out of specified interval.
/// </summary>
/// <param name="min - Minimum value to be possible random generated"></param>
/// <param name="max - Maximum value to be possible random generated"></param>
/// <returns>Maximum diagonal value</returns>
float Matrix::Generate(double min, double max)
{
    if( min > max ) {
        double temp = max;
        max = min;
        min = max;
    }

    float fSumAbs = 0.0, fDiagElemAbs = 0.0, fDiff = 0.0, fMaxAii = FLT_MIN;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    uniform_real_distribution<double> dDistr(min, max);
    uniform_int_distribution<int> iDistr(1, 10);
    for( int i = 0; i < _iR; i++ )
    {
        fSumAbs = 0;
        for( int j = 0; j < _iC; j++ )
        {
            _p[i][j] = dDistr(gen);
            fSumAbs += abs(_p[i][j]);
            if( i == j )
            {
                fDiagElemAbs = abs(_p[i][j]);
                fSumAbs -= fDiagElemAbs;
            }
        }
        if( fDiagElemAbs <= fSumAbs )
        {
            // if we don't have dominance of the diagonal
            fDiff = fSumAbs - fDiagElemAbs;
            if( _p[i][i] >= 0 )
            {
                _p[i][i] += (fDiff + iDistr(gen));
            }
            else
            {
                _p[i][i] -= (fDiff + iDistr(gen));
            }
        }

    }
    for( int i = 0; i < _iR; i++ ) {
        if( _p[i][i] > fMaxAii )
            fMaxAii = _p[i][i];
    }
    return fMaxAii;
}
/// <summary>
/// Print matrix to comand line
/// </summary>
void Matrix::PrintMatrixToShell()
{
    int width = 12;
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

/*
//----------------------------------------------------------------------------------------

	//int     rows = 3; int cols = rows;
	//// Init all necessary vectors and matrix A
	//Matrix  A = Matrix(rows, cols);
	//float   fMaxDiagVal = A.Generate(-50, 50);
	//Vector  vXZero = Vector(rows);
	//Vector  b = Vector(rows);

	//b.Generate(-10.0, 10.0);
	//// calculate s_prmCzbsz - step 0
	//s_prmCzbsz.fBeta        = 2 * fMaxDiagVal;
	//s_prmCzbsz.fOmegaZero   = (s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha) / (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha);
	//s_prmCzbsz.fC           = 2 / (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha);
	//if( fabs(s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha) > 0.00001f )
	//{
	//    s_prmCzbsz.fL = 2 * (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha) / (s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha);
	//}
	//else {
	//    s_prmCzbsz.fL = FLT_MAX;
	//    printf ("L value is +inf!!!");
	//}
	// end of step 0

	//TEST
	/*
	A.PrintMatrixToShell();
	printf("Maximum diagonal value: %.6f\n", fMaxDiagVal);
	printf("Vector b:\n");
	b.PrintVectorToShell();
	printf("Vector x:\n");
	vXPrev.PrintVectorToShell();
	*/
	// END TEST
