#include "TestData.h"
#include "Equation.h"
#include "Matrix.h"

TestData::TestData()
{
	this->data_file_paths[10] = "equations/equation_data_10.txt";
	this->data_file_paths[30] = "equations/equation_data_30.txt";
	this->data_file_paths[50] = "equations/equation_data_50.txt";
	this->data_file_paths[70] = "equations/equation_data_70.txt";
	this->data_file_paths[100] = "equations/equation_data_100.txt";
	this->data_file_paths[150] = "equations/equation_data_150.txt";
	this->data_file_paths[200] = "equations/equation_data_200.txt";
	this->data_file_paths[300] = "equations/equation_data_300.txt";
	this->data_file_paths[400] = "equations/equation_data_400.txt";
	this->data_file_paths[500] = "equations/equation_data_500.txt";
	this->data_file_paths[600] = "equations/equation_data_600.txt";
	this->data_file_paths[700] = "equations/equation_data_700.txt";
}

std::string TestData::getDataFilePath(int matrix_size)
{
	return this->data_file_paths[matrix_size];
}

void TestData::generateEquation(int matrix_size, int min_val, int max_val, std::string file_path)
{
	// Init vactor b and matrix A
	int rows = matrix_size;
	int cols = rows;

	Matrix  A = Matrix(rows, cols);
	Vector  b = Vector(rows);

	// generate all necessary data 
	float   fMaxDiagVal = A.Generate(min_val, max_val);
	b.Generate(min_val, max_val);

	// save generated data to file 
	Equation eq = Equation();
	eq.PrintEquationToFile(A, b, file_path);
}

void TestData::generateDataForTests()
{
	int matrix_size = 0;

	matrix_size = 10;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 30;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 50;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 70;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 100;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 150;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 200;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 300;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 400;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 500;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 600;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
	matrix_size = 700;
	generateEquation(matrix_size, -100, 100, getDataFilePath(matrix_size));
}
