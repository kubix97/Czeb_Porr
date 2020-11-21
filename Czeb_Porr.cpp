// Czeb_Porr.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <random>
#include <chrono>

using namespace std;

class Matrix 
{
    int _r, _c;
    float** _p;
    
    public:
        /**
         * Allocate memory block for our matrix, without initialization
         @param rows Number of rows in matrix
         @param cols Number of columns in matrix
        */
        Matrix(int rows, int cols) 
        {
            _r = rows; _c = cols;

            /*_p = new float*[_r];
            float *temp = new float [_r * _c];
            for( int i = 0; i < _r; i++ ) {
                _p[i] = temp + i * _c;
            }*/

            _p = (float **)malloc( _r * sizeof(float*) + _r *_c*sizeof(float) );
            float *temp = (float*) (_p + _r);
            for( int i = 0; i < _r; i++ ) {
                _p[i] = temp + i * _c;
            }
        }

        ~Matrix() 
        {
            free( _p );
        }

        int GetNumOfCals() { return _c; }
        int GetNumOfRows() { return _r; }
        
        /**
         * Fills matrix with random float values from interval specified by params. 
         * Ensure the dominance of diagonal, 
           to meet this condition values from diagonal 
           can be out of specified interval.
         * @param min Minimum value to be possible random generated
         * @param max Maximum value to be possible random generated
         * @return Maximum diagonal value
        */
        float Generate(double min, double max) 
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
            for( int i = 0; i < _r; i++ ) 
            {
                fSumAbs = 0;
                for( int j = 0; j < _c; j++ ) 
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
            for( int i = 0; i < _r; i++ ) {
                if( _p[i][i] > fMaxAii )
                    fMaxAii = _p[i][i];
            }
            return fMaxAii;
        }

        void PrintMatrixToShell()
        {
            int width = 12;
            printf("Matrix:\n");
            for( int i = 0; i < _r; i++ )
            {
                for( int j= 0; j < _c-1; j++ ) 
                {
                    printf( "%*.6f", width, _p[i][j] );
                }
                printf("%*.6f\n", width, _p[i][_c-1]);
            }
        }

        float* operator [] (int i) { return _p[i]; }
};

class Vector
{
    int _r;
    float* _pv;

    public:
        /**
         * Initialize column vector filled with zero values
         @param lenght
        */
        Vector( int length) 
        {
            _r = length;
            _pv = (float*) malloc( _r * sizeof(float) );
            for( int i = 0; i < _r; i++ ) 
                _pv[i] = 0.0;

            // We can use the calloc instead,
            // but due to later works in project we decide to use malloc
            //_pv = (float*) calloc(_r, sizeof(float));
        }

        ~Vector()
        {
            free(_pv);
        }

        void PrintVectorToShell()
        {
            for( int i = 0; i < _r-1 ; i++ )
            {
                printf( "%.6f\t", _pv[i] );
            }
            printf("%.6f\n", _pv[_r-1]);
        }

        /**
         Filling vector with random float values
        */
        void Generate(double min, double max)
        {
            if( min > max ) {
                double temp = max;
                max = min;
                min = max;
            }
            unsigned seed = chrono::system_clock::now().time_since_epoch().count();
            default_random_engine gen(seed);
            uniform_real_distribution<double> dDistr(min, max);
            for( int i = 0; i < _r; i++ )
                _pv[i] = dDistr(gen);

        }

};

int main()
{
    int rows = 3; int cols = rows;
    Matrix A = Matrix(rows,cols);
    float max = A.Generate(-50, 50);
    A.PrintMatrixToShell();
    printf("Maximum diagonal value: %.6f\n", max);
    Vector x = Vector(rows);
    printf("Vector x:\n");
    x.PrintVectorToShell();
    Vector b = Vector(rows);
    printf("Vector b:\n");
    b.Generate(-10.0, 10.0);
    b.PrintVectorToShell();
}