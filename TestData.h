#ifndef TESTDATA_H
#define TESTDATA_H

#include <iostream>
#include <map>

class TestData
{
private:
	std::map<int, std::string> data_file_paths;	// przechowuj sciezki do plikow z macierza o wielkosci x jako wartosc dla klucza rownego wielkosci macierzy

	/*
	 * Metoda generujaca rownanie do testow.
	 *
	 * Metoda generuje macierz wspolczynnikow A oraz wektor b wyrazow wolnych oraz zapisuje
	 * wygenerowane dane do pliku.
	 *
	 * @param	int			matrix_size					- rozmiar macierzy jaka wygenerujemy do testu
	 * @param	int			min_val						- minimalna wartosc wspolczynnikow w macierzy
	 * @param	int			max_val						- maksymalna wartosc wspolczynnikow w macierzy
	 * @param	std::string file_path					- sciezka do pliku, w ktorym zapiszemy dane rownania
	 *
	 */
	void generateEquation(int matrix_size, int min_val, int max_val, std::string file_path);
public:

	/*
	 * Konstruktor. 
	 * W konstruktorze wypelniamy mape >data_file_paths
	 */
	TestData();

	/* 
	 * Zwroc sciezke do pliku dla klucza bedacego rozmiarem macierzy zapisanej w danym pliku 
	 */
	std::string getDataFilePath(int matrix_size);

	/*
	 * Metoda generujaca rownania potrzebne do przeprowadzenia testow.
	 */
	void generateDataForTests();
};

#endif // !TESTDATA_H