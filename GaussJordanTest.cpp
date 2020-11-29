#include "GaussJordanTest.h"
#include "GaussJordan.h"
#include "Matrix.h"
#include "Vector.h"
#include "Equation.h"

void GaussJordanTest::writeTestResultsToFile(int matrix_size, double results[], std::string result_file_path)
{
	std::ofstream file;
	file.open(result_file_path, std::ios_base::app);
	if (file.is_open())
	{
		//file << "#Gauss-Jordan elimination execution time test" << "\n\n";
		/*file << std::fixed << std::setprecision(9);
		file << "Results for " << matrix_size << " rows" << "\n";
		file << "\t" << "Wersja sekwencyjna:"		<<	"\n\t\t"	<<	results[0] * 1e3	<< " ms\n";
		file << "\t" << "Wersja zrównoleglona:"		<<	"\n\t\t"	<<	results[1] * 1e3	<< " ms\n";
		file << "\t" << "Wersja z wektoryzacj¹:"	<<	"\n\t\t"	<<	results[2] * 1e3	<< " ms\n";*/

		file << std::fixed << std::setprecision(9);
		file << matrix_size << " rows:" << "\n";
		file << "\t" << "Wersja sekwencyjna:" << "\t\t\t" << results[0] * 1e3 << " ms\n";
		file << "\t" << "Wersja zrównoleglona:" << "\t\t" << results[1] * 1e3 << " ms\n";
		file << "\t" << "Wersja z wektoryzacj¹:" << "\t\t" << results[2] * 1e3 << " ms\n\n";
		
		file.close();
	}
	else {
		printf("Unable to wtrite to file");
	}
}

GaussJordanTest::GaussJordanTest()
{
	this->result_file_path = "GaussJordanEliminationTestResult.txt";
	this->temp_equation_file_path = "GaussJordanEliminationTestTempEq.txt";
}

GaussJordanTest::GaussJordanTest(std::string result_file, std::string temp_equation_file)
{
	this->result_file_path = result_file;
	this->temp_equation_file_path = temp_equation_file;
}

void GaussJordanTest::performTest(int matrix_size, int min_val, int max_val, std::string file_path, std::string result_file_path, int number_of_executions)
{
	// Init all necessary vectors and matrix A
	double results[3] = {};
	int rows = matrix_size;
	int cols = rows;

	Matrix  A = Matrix(rows, cols);
	Vector  vXZero1 = Vector(rows);
	Vector  vXZero2 = Vector(rows);
	Vector  vXZero3 = Vector(rows);
	Vector  b = Vector(rows);

	// generate all necessary data 
	float   fMaxDiagVal = A.Generate(min_val, max_val);
	b.Generate(min_val, max_val);

	// save generated data to file 
	Equation eq = Equation();
	eq.PrintEquationToFile(A, b, file_path);

	// init necessary matrix A and vector b for each kind of calculation
	Matrix A1 = Matrix(rows, cols);
	Vector b1 = Vector(rows);
	Matrix A2 = Matrix(rows, cols);
	Vector b2 = Vector(rows);
	Matrix A3 = Matrix(rows, cols);
	Vector b3 = Vector(rows);

	for (int i = 0; i < number_of_executions; i++) {
		// read data from file
		eq.ReadEquationFromFile(A1, b1, file_path);
		eq.ReadEquationFromFile(A2, b2, file_path);
		eq.ReadEquationFromFile(A3, b3, file_path);


		// perform calculation
		std::cout << "\n\nExecution time for " << rows << " rows" << "\n\n";
		std::cout << "Sekwencyjna:\t\t " << std::endl;
		results[0] = GaussJordan::GaussJordanAlgotithm(A1, vXZero1, b1);
		std::cout << "Zrownoleglona:\t\t " << std::endl;
		results[1] = GaussJordan::GaussJordanAlgorithmWithParalelization(A2, vXZero2, b2);
		std::cout << "Wektoryzacja:\t\t " << std::endl;
		results[2] = GaussJordan::GaussJordanAlgorithmWithVectorization(A3, vXZero3, b3);

		// write execution time of each method to file
		writeTestResultsToFile(rows, results, result_file_path);
	}
}

void GaussJordanTest::performTestGroup(int test_repetition_number)
{
	GaussJordanTest test = GaussJordanTest();
	test.performTest(15,	-50,	50,	"gaussjordantest.txt",	this->result_file_path, 5);
	test.performTest(40,	-50,	50,	"gaussjordantest.txt",	this->result_file_path, 5);
	test.performTest(100,	-50,	50,	"gaussjordantest.txt",	this->result_file_path, 5);
	test.performTest(300,	-50,	50,	"gaussjordantest.txt",	this->result_file_path, 5);
	test.performTest(500,	-50,	50,	"gaussjordantest.txt",	this->result_file_path, 5);
}
