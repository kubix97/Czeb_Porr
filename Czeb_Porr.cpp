// Czeb_Porr.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <random>
#include <chrono>
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


Vector CzebAlg(Matrix &A, Vector &x0, Vector &b, int iters, int s)
{
    Vector  vXPrev = Vector(x0.GetLen());
    int     k = 0;
    Vector  vX = Vector(x0.GetLen()); // current x start filled with zeros
    Vector  vTemp, vTm;
    float   fCff;

    for( int i = 0; i < iters; i++ )
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
    }
    return vX;
}

int main()
{
    int     rows = 3; int cols = rows;
    // Init all necessary vectors and matrix A
    Matrix  A = Matrix(rows, cols);
    float   fMaxDiagVal = A.Generate(-50, 50);
    Vector  vXZero = Vector(rows);
    Vector  b = Vector(rows);

    b.Generate(-10.0, 10.0);
    // calculate s_prmCzbsz - step 0
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


    
    
    
}