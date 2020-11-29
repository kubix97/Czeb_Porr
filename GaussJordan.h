#ifndef GAUSS_H
#define GAUSS_H

#include <chrono>
#include "Matrix.h"
#include "Vector.h"
class GaussJordan
{
public:
	/**
	  * Statyczna metoda obliczajaca rozwiazania ukladu rownan liniowych 
	  * zdefioniowanego wzorem Ax = b, za pomoca metody eliminacji Gaussa-Jordana
	  *
	  * @param	Matrix &A	- macierz wspolczynnikow ukladu rownan
	  * @param	Vector &x	- wektor rozwiazan ukladu rownan
	  * @param	Vector &b	- wektor wyrazow wolnych ukladu rownan
	  * 
	  * @retval double		- czas wykonania obliczen w sekundach
	  */
	static double GaussJordanAlgotithm(Matrix &A, Vector &x , Vector &b);

	/**
	  * Metoda rozszerzajaca metode GaussJordanAlgotithm() o dyrektywe 
	  * zrownoleglajaca obliczenia
	  *
	  * @see GaussJordanAlgotithm()
	  * @param	Matrix &A	- macierz wspolczynnikow ukladu rownan
	  * @param	Vector &x	- wektor rozwiazan ukladu rownan
	  * @param	Vector &b	- wektor wyrazow wolnych ukladu rownan
	  *
	  * @retval double		- czas wykonania obliczen w sekundach
	  */
	static double GaussJordanAlgorithmWithParalelization(Matrix& A, Vector& x, Vector& b);

	/**
	  * Metoda rozszerzajaca metode GaussJordanAlgotithm() o dyrektywe
	  * wektoryzacji obliczen
	  *
	  * @see GaussJordanAlgotithm()
	  * @param	Matrix &A	- macierz wspolczynnikow ukladu rownan
	  * @param	Vector &x	- wektor rozwiazan ukladu rownan
	  * @param	Vector &b	- wektor wyrazow wolnych ukladu rownan
	  *
	  * @retval double		- czas wykonania obliczen w sekundach
	  */
	static double GaussJordanAlgorithmWithVectorization(Matrix& A, Vector& x, Vector& b);

	/**
	  * Statyczna metoda wypisujaca w konsoli wynik obliczen w postaci 
	  * wspolczynnikow wektora rozwiazan x
	  *
	  * @param	Vector &x	- wektor rozwiazan ukladu rownan
	  *
	  */
	static void PrintResults(Vector& x);
};

#endif // GAUSS_H