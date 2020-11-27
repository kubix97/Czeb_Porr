#include <iostream>
#include "Gauss.h"
#include "Matrix.h"
#include "Vector.h"

void Gauss::GaussJordanAlgotithm(Matrix &A, Vector &x, Vector &b)
{
	float alfa;
	int i = 0, j = 0, k = 0;
	int n = A.GetNumOfRows() - 1;					// -1, because indexing starts from 0
	std::cout << "Wymiar macierzy wynosi: "
		<< A.GetNumOfRows() << std::endl;
	std::cout.precision(6);								// set precision
	std::cout.setf(std::ios::fixed);

	A.PrintMatrixToShell();

	/* calculate Gauss-Jordan elimination */
	for (k = 0; k <= n; k++) {						// petla wierszy eliminujacych
													// (kolumn eliminowanych)

		for (i = 0; i <= n; i++) {					// petla modyfikacji wierszy
													// ponizej i powyzej
			if (i == k) continue;
			alfa = A[i][k] / A[k][k];
			for (j = k; j <= n; j++)				// petla kolumn w danym wierszu
				A[i][j] = A[i][j] - alfa * A[k][j];
			b[i] = b[i] - alfa * b[k];
		}
	}

	A.PrintMatrixToShell();							// print matrix after Gauss-Jordan elimination	

	/* calculate x vector which contains results */
	for (k = 0; k <= n; k++)
		x[k] = b[k] / A[k][k];

	printf("\nResults:\n");
	for (i = 0; i <= n; i++)						// print results
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
}

void Gauss::GaussJordanAlgorithmWithParalelization(Matrix& A, Vector& x, Vector& b)
{
	float alfa;
	int i = 0, j = 0, k = 0;
	int n = A.GetNumOfRows() - 1;					// -1, because indexing starts from 0
	std::cout << "Wymiar macierzy wynosi: "
		<< A.GetNumOfRows() << std::endl;
	std::cout.precision(6);								// set precision
	std::cout.setf(std::ios::fixed);

	A.PrintMatrixToShell();

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

	A.PrintMatrixToShell();							// print matrix after Gauss-Jordan elimination	

	/* calculate x vector which contains results */
	for (k = 0; k <= n; k++)
		x[k] = b[k] / A[k][k];

	printf("\nResults:\n");
	for (i = 0; i <= n; i++)						// print results
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
}

void Gauss::GaussJordanAlgorithmWithVectorization(Matrix& A, Vector& x, Vector& b)
{
	float alfa;
	int i = 0, j = 0, k = 0;
	int n = A.GetNumOfRows() - 1;					// -1, because indexing starts from 0
	std::cout << "Wymiar macierzy wynosi: "
		<< A.GetNumOfRows() << std::endl;
	std::cout.precision(6);								// set precision
	std::cout.setf(std::ios::fixed);

	A.PrintMatrixToShell();

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

	A.PrintMatrixToShell();							// print matrix after Gauss-Jordan elimination	

	/* calculate x vector which contains results */
	for (k = 0; k <= n; k++)
		x[k] = b[k] / A[k][k];

	printf("\nResults:\n");
	for (i = 0; i <= n; i++)						// print results
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
}
