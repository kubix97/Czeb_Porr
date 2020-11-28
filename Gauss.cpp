#include <iostream>
#include <chrono>
#include <omp.h>
#include "Gauss.h"
#include "Matrix.h"
#include "Vector.h"

void Gauss::GaussJordanAlgotithm(Matrix &A, Vector &x, Vector &b)
{
	float alfa;
	int i = 0, j = 0, k = 0;
	int n = A.GetNumOfRows() - 1;								// -1, because indexing starts from 0

	//A.PrintMatrixToShell();

	auto start = std::chrono::high_resolution_clock::now();		// start time measurment

	/* calculate Gauss-Jordan elimination */
	for (k = 0; k <= n; k++) {									// petla wierszy eliminujacych
																// (kolumn eliminowanych)

		for (i = 0; i <= n; i++) {								// petla modyfikacji wierszy
																// ponizej i powyzej
			if (i == k) continue;
			alfa = A[i][k] / A[k][k];
			for (j = k; j <= n; j++)							// petla kolumn w danym wierszu
				A[i][j] = A[i][j] - alfa * A[k][j];
			b[i] = b[i] - alfa * b[k];
		}
	}

	//A.PrintMatrixToShell();									// print matrix after Gauss-Jordan elimination	

	/* calculate x vector which contains results */
	for (k = 0; k <= n; k++)
		x[k] = b[k] / A[k][k];

	auto end = std::chrono::high_resolution_clock::now();		// stop time measurment
	double duration = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	duration *= 1e-9;
	printf("Execution time = %.9f sec = %.9f ms\n", duration, duration * 1e3);
}

void Gauss::GaussJordanAlgorithmWithParalelization(Matrix& A, Vector& x, Vector& b)
{
	float alfa;
	int i = 0, j = 0, k = 0;
	int n = A.GetNumOfRows() - 1;					// -1, because indexing starts from 0

	//A.PrintMatrixToShell();

	auto start = std::chrono::high_resolution_clock::now();

	/* calculate Gauss-Jordan elimination */
	for (k = 0; k <= n; k++) {						// petla wierszy eliminujacych
													// (kolumn eliminowanych)

		#pragma omp parallel for private(alfa,j)	// dyrektywa zrownoleglajaca
		for (i = 0; i <= n; i++) {					// petla modyfikacji wierszy
													// ponizej i powyzej
			if (i == k) continue;
			alfa = A[i][k] / A[k][k];
			for (j = k; j <= n; j++)				// petla kolumn w danym wierszu
				A[i][j] = A[i][j] - alfa * A[k][j];
			b[i] = b[i] - alfa * b[k];
		}
	}

	//A.PrintMatrixToShell();							// print matrix after Gauss-Jordan elimination	

	/* calculate x vector which contains results */
	for (k = 0; k <= n; k++)
		x[k] = b[k] / A[k][k];

	auto end = std::chrono::high_resolution_clock::now();
	double duration = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	duration *= 1e-9;
	printf("Execution time = %.9f sec = %.9f ms\n", duration, duration * 1e3);
}

void Gauss::GaussJordanAlgorithmWithVectorization(Matrix& A, Vector& x, Vector& b)
{
	float alfa;
	int i = 0, j = 0, k = 0;
	int n = A.GetNumOfRows() - 1;					// -1, because indexing starts from 0

	//A.PrintMatrixToShell();

	auto start = std::chrono::high_resolution_clock::now();

	/* calculate Gauss-Jordan elimination */
	for (k = 0; k <= n; k++) {						// petla wierszy eliminujacych
													// (kolumn eliminowanych)

		#pragma omp simd							// dyrektywa zrownoleglajaca
		for (i = 0; i <= n; i++) {					// petla modyfikacji wierszy
													// ponizej i powyzej
			if (i == k) continue;
			alfa = A[i][k] / A[k][k];
			for (j = k; j <= n; j++)				// petla kolumn w danym wierszu
				A[i][j] = A[i][j] - alfa * A[k][j];
			b[i] = b[i] - alfa * b[k];
		}
	}

	//A.PrintMatrixToShell();							// print matrix after Gauss-Jordan elimination	

	/* calculate x vector which contains results */
	for (k = 0; k <= n; k++)
		x[k] = b[k] / A[k][k];

	auto end = std::chrono::high_resolution_clock::now();
	double duration = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	duration *= 1e-9;
	printf("Execution time = %.9f sec = %.9f ms\n", duration, duration * 1e3);
}

void Gauss::PrintResults(Vector& x)
{
	printf("\nResults:\n");
	for (int i = 0; i < x.GetLen(); i++)
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
}
