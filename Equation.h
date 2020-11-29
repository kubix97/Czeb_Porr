#ifndef EQUATION_H
#define EQUATION_H

#include <fstream>
#include <iomanip>
#include "Matrix.h"
#include "Vector.h"

class Equation
{
public:
	Equation() {}
	~Equation() {}

	/*
	 * Metoda zapisujaca do pliku wygenerowana macierz wspolczynnikow 
	 * ukladu rownan - A oraz wektor wyrazow wolnych - b
	 *
	 * @param	Matrix &A				-	macierz wspolczynnikow ukladu rownan
	 * @param	Vector &b				-	wektor wyrazow wolnych ukladu rownan
	 * @param	std::string file_path	-	sciezka do pliku, w ktorym maja zostac zapisane wygenerowane dane,
	 *										domyslnie - "Dane.txt"
	 * 
	 */
	void PrintEquationToFile(Matrix& A, Vector& b, std::string file_path);

	void PrintEquationToFile(Matrix& A, Vector& b);

	
	/*
	 * Metoda wczytujaca z pliku macierz wspolczynnikow ukladu rownan - A
	 * oraz wektor wyrazow wolnych - b
	 *
	 * @param	Matrix		&A			-	macierz wspolczynnikow ukladu rownan
	 * @param	Vector		&b			-	wektor wyrazow wolnych ukladu rownan
	 * @param	std::string file_path	-	sciezka do pliku, z ktorym maja zostac wczytane dane,
	 *										domyslnie - "Dane.txt"
	 *
	 */
	void ReadEquationFromFile(Matrix& A, Vector& b, std::string file_path);

	void ReadEquationFromFile(Matrix& A, Vector& b);

};

#endif