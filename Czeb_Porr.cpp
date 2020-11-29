// Czeb_Porr.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <iomanip>
#include "Matrix.h"
#include "Gauss.h"


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

void PrintEquationToFile(Matrix& A, Vector& b)
{
	ofstream file;
	int iR = A.GetNumOfRows(), iC = A.GetNumOfCals();
	file.open("Dane.txt");
	if (file.is_open())
	{
		file << fixed << setprecision(6);
		file << iR << "\t" << iC << "\n";
		for (int i = 0; i < iR; i++)
		{
			for (int j = 0; j < iC; j++) {
				file << A[i][j] << "\t";
			}
			file << "\n";
		}
		for (int i = 0; i < iR; i++) {
			file << b[i] << "\t";
		}
		file.close();
	}
	else {
		printf("Unable to wtrite to file");
	}
}

void ReadEquationFromFile(Matrix& A, Vector &b)
{
	ifstream file;
	int iR = 0, iC = 0;

	file.open("Dane.txt");
	if (file.is_open())
	{
		file >> fixed >> setprecision(6);
		file >> iR >> iC;

		// check if matrix has proper size
		if (iR <= 0 && iC <= 0) {
			printf("Matrix must have at least one row");
			return;
		}

		// read A matrix
		for (int i = 0; i < iR; i++)
			for (int j = 0; j < iC; j++)
				file >> A[i][j];
		
		// read b vector
		for (int i = 0; i < iC; i++)
			file >> b[i];

		file.close();
	}
	else {
		printf("Unable to wtrite to file");
	}
}

int main()
{
	int rows = 700; 
	int cols = rows;

	// Init all necessary vectors and matrix A
	Matrix  A = Matrix(rows, cols);
	float   fMaxDiagVal = A.Generate(-50, 50);
	Vector  vXZero1 = Vector(rows);
	Vector  vXZero2 = Vector(rows);
	Vector  vXZero3 = Vector(rows);
	Vector  b = Vector(rows);
	b.Generate(-50, 50);

	PrintEquationToFile(A, b);
	

	Matrix A1 = Matrix(rows, cols);
	Vector b1 = Vector(rows);
	ReadEquationFromFile(A1, b1);
	/*A1.PrintMatrixToShell();
	cout << "Wektor b1: " << endl;
	b1.PrintVectorToShell();*/

	Matrix A2 = Matrix(rows, cols);
	Vector b2 = Vector(rows);
	ReadEquationFromFile(A2, b2);
	/*A2.PrintMatrixToShell();
	cout << "Wektor b2: " << endl;
	b2.PrintVectorToShell();*/

	Matrix A3 = Matrix(rows, cols);
	Vector b3 = Vector(rows);
	ReadEquationFromFile(A3, b3);
	/*A3.PrintMatrixToShell();
	cout << "Wektor b3: " << endl;
	b3.PrintVectorToShell();*/

	
	//Matrix A1 = A;
	//b.PrintVectorToShell();
	//Vector b1(b);
	//b1.PrintVectorToShell();

	//Matrix A2 = A;
	//Vector b2(b);
	//b2.PrintVectorToShell();

	cout << "Execution time for " << rows << " rows" << "\n\n";

	cout << "Sekwencyjna:\t\t " << endl;
	Gauss::GaussJordanAlgotithm(A1, vXZero1, b1);
	cout << "Zrownoleglona:\t\t " << endl;
	Gauss::GaussJordanAlgorithmWithParalelization(A2, vXZero2, b2);
	cout << "Wektoryzacja:\t\t " << endl;
	Gauss::GaussJordanAlgorithmWithVectorization(A3, vXZero3, b3);

	/*cout << "Sekwencyjna:\t\t " << endl;
	Gauss::PrintResults(vXZero1);
	cout << "Zrownoleglona:\t\t " << endl; 
	Gauss::PrintResults(vXZero2);
	cout << "Wektoryzacja:\t\t " << endl;
	Gauss::PrintResults(vXZero3);
	A3.PrintMatrixToShell();*/

	return 0;
    
}