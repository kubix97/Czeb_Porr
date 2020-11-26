// Czeb_Porr.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <random>
#include <chrono>
#include <iomanip>
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
	int n = 2;									// rozmiar macierzy
	int i = 0, j = 0, k = 0;					// zmienne pomocnicze
	float alfa;
	cout.precision(4);							// ustalenie precyzji
	cout.setf(ios::fixed);
	//float a[3][3] = { {1,  1, -1},
	//				  {1, -1,  2},
	//				  {2,  1,  1} };
	//float b[3] = { 7,3,9 };
	float x[3] = {};
	float a[3][3] = {	{6,   2, -1},
						{3,   9,  1},
						{2,  -2,  3} };
	float b[3] = { 5,9,-1 };

	for (i = 0; i <= n; i++)					//print the new matrix
	{
		for (j = 0; j <= n; j++)
			cout << a[i][j] << setw(16);
		cout << "\n";
	}
	for (k = 0; k <= n; k++) {					// petla wierszy eliminujacych
												// (kolumn eliminowanych)
		#pragma omp parallel for private(alfa,j)
		for (i = 0; i <= n; i++) {				// petla modyfikacji wierszy
												// ponizej i powyzej
			if (i == k) continue;
			alfa = a[i][k] / a[k][k];
			for (j = k; j <= n; j++)			// petla kolumn w danym wierszu
				a[i][j] = a[i][j] - alfa * a[k][j];
			b[i] = b[i] - alfa * b[k];
		}
	}

	cout << "\n";
	for (i = 0; i <= n; i++)					// macierz po eliminacji Gaussa - Jordana
	{
		for (j = 0; j <= n; j++)
			cout << a[i][j] << setw(16);
		cout << "\n";
	}

	for (k = 0; k <= n; k++) 
		x[k] = b[k] / a[k][k];


	cout << "\nThe values of result vector are:\n";
	for (i = 0; i <= n; i++)
		cout << "x[" << i << "] = " << x[i] << endl;				// Wypisz rozwiazanie  


	return 0;
    
}