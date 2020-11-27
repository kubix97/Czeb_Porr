// Czeb_Porr.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include "Matrix.h"


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


Vector CzebAlg(Matrix &A, Vector &x0, Vector &b, int maxIters, int s, float accuracy)
{
    Vector  vXPrev = Vector(x0.GetLen());
    int     k = 0;
    Vector  vX = Vector(x0.GetLen()); // current x start filled with zeros
    Vector  vTemp, vTm;
    float   fCff;
    auto start = chrono::high_resolution_clock::now();
    for( int i = 0; i < maxIters; i++ )
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
            printf("Cheb iteration stops after %d iterations\n", i);
            maxIters = i; // save num of iterations to calcuate execution time per iteration
            break;
        }
            
    }
    auto end = chrono::high_resolution_clock::now();
    double duration = (double) chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    duration *= 1e-9;
    printf("Execution time = %.9f sec = %.9f ms\nAverage time per iteration = %f.9 ms\n", duration, duration*1e3, duration*1e3/maxIters);
    return vX;
}

void PrintEquationToFile(Matrix& A, Vector& b)
{
    ofstream file;
    int iR = A.GetNumOfRows(), iC = A.GetNumOfCals();
    file.open("Dane.txt");
    if( file.is_open() )
    {
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

int main()
{
    //int     rows = 3; int cols = rows;
    int     rows = 5; int cols = rows;

    // TEST
    
    /*float tabAtst[5][5] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 1.0f};*/
    float tabAtst[5][5] = {1.0f, 2.0f, 3.0f, 4.0f, 1.0f,
        -1.0f, -4.0f, 3.0f, 4.0f, 1.0f,
        2.0f, 3.0f, 4.0f, -5.0f, 1.0f,
        -6.0f, 3.0f, 2.0f, 1.0f, -2.0f,
        6.0f, -2.0f, -1.0f, 1.0f, 4.0f};

    Matrix Atst = Matrix(rows,cols);
    Atst.GenWithFixedVal(tabAtst);
    Atst.PrintMatrixToShell();
    Matrix Test= Atst.Transpose();
    //Test.PrintMatrixToShell();
    Atst = Atst.MultiplyWithMatrix(Test);
    Atst.PrintMatrixToShell();
    //float tabB[4] = {2.0f, 2.0f, 2.0f, 2.0f};
    float tabB[5] = {-5.0f, 51.0f, 17.0f, -165.0f, 168.0f};
    //float tabC[5] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    Vector btst = Vector(rows);
    btst.GenWithFixedVal(rows, tabB);
    printf("Vector bTest:\n");
    btst.PrintVectorToShell();
    PrintEquationToFile(Atst, btst);
    float fMaxDiagVal = Atst.GetMaxDiagVal();
    printf("Max diagonal value: %f\n", fMaxDiagVal);

   /* vSub.GenWithFixedVal(rows, tabC);   
    printf("Vector vSub:\n");
    vSub.PrintVectorToShell();
    Vector vRes = Atst.MultiplyWithVector(btst);
    printf("Vector vRes:\n");
    vRes.PrintVectorToShell();
    vRes = vRes.Add(vSub);
    printf("Vector vResMod:\n");
    vRes.PrintVectorToShell();
    vRes = vRes.MultiplyByVal(2);
    printf("Vector vResMod:\n");
    vRes.PrintVectorToShell();*/

    Vector  vXZero = Vector(rows);
    // calculate s_prmCzbsz
    s_prmCzbsz.fBeta        = 2 * fMaxDiagVal;
    s_prmCzbsz.fOmegaZero   = (s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha) / (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha);
    s_prmCzbsz.fC           = 2 / (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha);
    if( fabs(s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha) > 0.00001f )
    {
        s_prmCzbsz.fL = 2 * (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha) / (s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha);
    }
    else {
        s_prmCzbsz.fL = FLT_MAX;
        printf ("L value is +inf!!!");
    }
    Vector Xtst = CzebAlg( Atst, vXZero, btst, 1000, 30000, 1e-8f);
    printf("Vector Xres:\n");
    Xtst.PrintVectorToShell(); // X result should be [-1;2;1;-2;1]
    // END TEST

    // Init all necessary vectors and matrix A
    //Matrix  vA = Matrix(rows, cols);
    //float   fMaxDiagVal = vA.Generate(-50, 50);
    //Vector  vXZero = Vector(rows);
    //Vector  vB = Vector(rows);
    //vB.Generate(-100.0, 100.0);
    //// calculate s_prmCzbsz
    //s_prmCzbsz.fBeta = 2 * fMaxDiagVal;
    //s_prmCzbsz.fOmegaZero = (s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha) / (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha);
    //s_prmCzbsz.fC = 2 / (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha);
    //if( fabs(s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha) > 0.00001f )
    //{
    //    s_prmCzbsz.fL = 2 * (s_prmCzbsz.fBeta + s_prmCzbsz.fAlpha) / (s_prmCzbsz.fBeta - s_prmCzbsz.fAlpha);
    //}
    //else {
    //    s_prmCzbsz.fL = FLT_MAX;
    //    printf("L value is +inf!!!");
    //}
    //Vector vXRes = CzebAlg(vA, vXZero, vB, 10000, 50000);
    //vXRes.PrintVectorToShell();
    //vA.PrintMatrixToShell();
    //vB.PrintVectorToShell();

    
    
    
}