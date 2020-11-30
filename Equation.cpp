#include "Equation.h"

void Equation::PrintEquationToFile(Matrix& A, Vector& b)
{
	std::ofstream file;
	int iR = A.GetNumOfRows(), iC = A.GetNumOfCals();
	file.open("Dane.txt");
	if (file.is_open())
	{
		file << std::fixed << std::setprecision(6);
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
		printf("Unable to write to file");
	}
}

void Equation::PrintEquationToFile(Matrix& A, Vector& b, std::string file_path)
{
	std::ofstream file;
	int iR = A.GetNumOfRows(), iC = A.GetNumOfCals();
	file.open(file_path);
	if (file.is_open())
	{
		file << std::fixed << std::setprecision(6);
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

void Equation::ReadEquationFromFile(Matrix& A, Vector& b)
{
	std::ifstream file;
	int iR = 0, iC = 0;

	file.open("Dane.txt");
	if (file.is_open())
	{
		file >> std::fixed >> std::setprecision(6);
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

void Equation::ReadEquationFromFile(Matrix& A, Vector& b, std::string file_path)
{
	std::ifstream file;
	int iR = 0, iC = 0;

	file.open(file_path);
	if (file.is_open())
	{
		file >> std::fixed >> std::setprecision(6);
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
