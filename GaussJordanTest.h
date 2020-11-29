#ifndef GAUSSJORDANTEST_H
#define GAUSSJORDANTEST_H

#include <iostream>

class GaussJordanTest
{
private:
	std::string result_file_path;
	std::string temp_equation_file_path;

	/*
	 * Metoda zapisujaca wynik testu do pliku
	 *
	 * @param	int				matrix_size			- rozmiar macierzy, dla ktorej ma zostac przeprowadzony test 
	 * @param	double			results[]			- tablica wynikow otrzymanych dla kazdej z metod obliczania
	 * @param	std::string		result_file_path	- sciezka do pliku, w ktorym maja zostac zapisane wyniki testu
	 *
	 */
	void writeTestResultsToFile(int matrix_size, double results[], std::string result_file_path);
public:
	GaussJordanTest();

	/*
	 * Konstruktor
	 *
	 * @param	std::string		result_file			-	sciezka do pliku, w ktorym maja zostac zapisane wyniki testow
	 * @param	std::string		temp_equation_file	-	sciezka do pliku, w ktorym maja byc przechowana wygenerowana 
	 *													macierz wspolczynnikow A oraz wektor wyrazow wolnych b
	 *
	 */
	GaussJordanTest(std::string result_file, std::string temp_equation_file);
	~GaussJordanTest() {}

	/*
	 * Metoda przeprowadzajaca test czasu wykonania eliminacji 
	 * Gaussa-Jordana.
	 * 
	 * Metoda generuje macierz wspolczynnikow A oraz wektor wyrazow wolnych.
	 *
	 * @param	int			matrix_size					- rozmiar macierzy jaka wygenerujemy do testu
	 * @param	int			min_val						- minimalna wartosc wspolczynnikow w macierzy
	 * @param	int			max_val						- maksymalna wartosc wspolczynnikow w macierzy
	 * @param	std::string file_path					- sciezka do pliku, w ktorym zapiszemy macierz
	 * @param	std::string result_file_path			- sciezka do pliku, w ktorym zapiszemy wyniki testu
	 * @param	int			number_of_executions		- liczba powtorzen testu na danej macierzy
	 *
	 */
	void performTest(int matrix_size, int min_val, int max_val, std::string file_path, std::string result_file_path, int number_of_executions);

	/*
	 * Metoda przeprowadzajaca scenariusz testowy dla zdefiniowanej wewnatrz grupy testow
	 *
	 * @param	int test_repetition_number	- liczba powtorzen kazdego ze zdefiniowanych w scenariuszu testow 
	 *
	 */
	void performTestGroup(int test_repetition_number);
};

#endif // !GAUSSJORDANTEST_H