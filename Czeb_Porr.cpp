// Czeb_Porr.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <iomanip>
#include "common.h"
#include "Matrix.h"

#include "GaussJordan.h"
#include "GaussJordanTest.h"
#include "Equation.h"
#include "TestData.h"

#if EXMODE == 2
#include <omp.h>
#endif


using namespace std;


struct CZBPRM {
    float fAlpha = 100.0;
    float fBeta = 0;
    float fOmegaZero = 0;
    float fOmega = 0;
    float fOmegaPrev = 0;
    float fC = 0;
    float fL = 0;
};

static CZBPRM s_prmCzbsz;

/// <summary>
/// 
/// </summary>
/// <param name="A"></param>
/// <param name="x0"></param>
/// <param name="b"></param>
/// <param name="maxIters"></param>
/// <param name="s"></param>
/// <param name="accuracy"></param>
/// <returns></returns>
Vector CzebAlg(Matrix &A, Vector &x0, Vector &b, int maxIters, int s, float accuracy)
{
    Vector  vXPrev = Vector(x0.GetLen());
    int     i, k = 0;
    Vector  vX = Vector(x0.GetLen()); // current x start filled with zeros
    Vector  vTemp, vTm;
    float   fCff;
    auto start = chrono::high_resolution_clock::now();
    for( i = 0; i < maxIters; i++ )
    {
        if( k >= s || i == 0 )
        {
            s_prmCzbsz.fOmega = s_prmCzbsz.fOmegaZero; s_prmCzbsz.fOmegaPrev = 0.0f; k = 0; vX = vXPrev;
        }
        vTemp = vX - vXPrev;
        vTemp.MultiplyByValI(s_prmCzbsz.fOmega * s_prmCzbsz.fOmegaPrev);
        vTemp.AddI(vX);

        fCff = s_prmCzbsz.fC * (1.0f + s_prmCzbsz.fOmega * s_prmCzbsz.fOmegaPrev);
        vTm  = A.MultiplyWithVector(vX) - b;
        vTm.MultiplyByValI(fCff);
        vTemp= vTemp - vTm; 

        vXPrev = vX; vX = vTemp; s_prmCzbsz.fOmegaPrev = s_prmCzbsz.fOmega; s_prmCzbsz.fOmega = 1 / (s_prmCzbsz.fL - s_prmCzbsz.fOmega);
        if( vX.CalcDistance(vXPrev) < accuracy ) {
            maxIters = i + 1; // save num of iterations to calcuate execution time per iteration
            break;
        }
            
    }
    auto end = chrono::high_resolution_clock::now();
    double duration = (double) chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    duration *= 1e-9;
    printf("\nCheb [MD: %d] iteration stops after %d iterations - matrix %d x %d\n", EXMODE, i + 1, A.GetNumOfRows(), A. GetNumOfCals());
    printf("Execution time = %.9f sec = %.9f ms\nAverage time per iteration = %.9f ms\n", duration, duration*1e3, duration*1e3/maxIters);
    return vX;
}

void PrintEquationToFile(Matrix& A, Vector& b)
{
    ofstream file;
    int iR = A.GetNumOfRows(), iC = A.GetNumOfCals();
    file.open("EquationsCheb//EquationCheb_10.txt");
    if( file.is_open() )
    {
        file << fixed << setprecision(6);
        file << iR << "\t" << iC << "\n";
        for( int i = 0; i < iR; i++ )
        {
            for( int j = 0; j < iC; j++ ) {
                file << A[i][j] << "\t";
            }
            file << "\n";
        }
        for( int i = 0; i < iR; i++ ) {
            file << b[i] << "\t";
        }
        file.close();
    }
    else { 
        printf("Unable to wtrite to file");
    }
}

void ReadEquationFromFile(Matrix& A, Vector& b)
{
    ifstream file;
    int iR = 0, iC = 0;

    file.open("equations//equation_data_10.txt");
    if( file.is_open() )
    {
        file >> fixed >> setprecision(6);
        file >> iR >> iC;

        // check if matrix has proper size
        if( iR <= 0 && iC <= 0 ) {
            printf("Matrix must have at least one row");
            return;
        }

        // read A matrix
        for( int i = 0; i < iR; i++ )
            for( int j = 0; j < iC; j++ )
                file >> A[i][j];

        // read b vector
        for( int i = 0; i < iC; i++ )
            file >> b[i];

        file.close();
    }
    else {
        printf("Unable to wtrite to file");
    }
}


#if EXMODE == 2
void OpenMpSet(int iDim) {
    //int iProcs = omp_get_num_procs();
    if( iDim <= 20 ) {
        //omp_set_dynamic(1);
        omp_set_num_threads(1);
    }
    else
    {
        omp_set_num_threads(4);
    }
}
#endif

/// <summary>
/// Czebyszew algorithm evaluation
/// </summary>
/// <returns></returns>
void Czebyszew(int iRows)
{
    int     rows = iRows; int cols = rows;
    bool    bReadFromFile = false;

#if EXMODE == 2
    OpenMpSet(rows);
#endif

    // Init all necessary vectors and matrix A

    Matrix A = Matrix(rows, cols);
    Vector vB = Vector(rows);
    if( !bReadFromFile )
    {
        Matrix  G = Matrix(rows, cols);
        G.Generate(-10, 10);
        Matrix Gt = G.Transpose();
        A = G.MultiplyWithMatrix(Gt);
        vB.Generate(-10000, 10000.0);
        //G.PrintMatrixToShell(); Gt.PrintMatrixToShell();
        PrintEquationToFile(A, vB);
    }
    else {
        ReadEquationFromFile(A, vB);
    }
    Vector  vXZero = Vector(rows);
    float fMaxDiagVal = A.GetMaxDiagVal();
    // calculate s_prmCzbsz
    s_prmCzbsz.fBeta = 2 * fMaxDiagVal;
    s_prmCzbsz.fOmegaZero = (s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha) / (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha);
    s_prmCzbsz.fC = 2 / (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha);
    if( fabs(s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha) > 0.00001f )
    {
        s_prmCzbsz.fL = 2 * (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha) / (s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha);
    }
    else {
        s_prmCzbsz.fL = FLT_MAX;
        printf("L value is +inf!!!");
    }
    Vector vXRes = CzebAlg(A, vXZero, vB, (int) 1e6, (int) 2e5, 6e-7f);
    //vXRes.PrintVectorToShell();
}

/// <summary>
/// Application starup
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(int argc, char* argv[])
{
    int     iMode = 0;
    if( argc >= 2 )
    {
        if( 0 == _stricmp(argv[1], "-gj") ) {
            iMode = 1;
        }
    }
    if( 0 == iMode )
    {
        Czebyszew(10);
        
        Czebyszew(30);
        

        Czebyszew(50);
        

        Czebyszew(70);
        

        Czebyszew(100);
       

        Czebyszew(150);
        

        Czebyszew(200);
        

        Czebyszew(300);
       

        Czebyszew(400);
       

        Czebyszew(500);
        

        Czebyszew(600);
        

        Czebyszew(700);
       

    }
    else
    {
        TestData td = TestData();
        GaussJordanTest test = GaussJordanTest();
        test.performTestGroupOnExistingEquation(5);
    }
}