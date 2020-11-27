#ifndef GAUSS_H
#define GAUSS_H

#include "Matrix.h"
#include "Vector.h"
class Gauss
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
	  */
	static void GaussJordanAlgotithm(Matrix &A, Vector &x , Vector &b);

	/**
	  * Metoda rozszerzajaca metode GaussJordanAlgotithm() o dyrektywe 
	  * zrownoleglajaca obliczenia
	  *
	  * @see GaussJordanAlgotithm()
	  * @param	Matrix &A	- macierz wspolczynnikow ukladu rownan
	  * @param	Vector &x	- wektor rozwiazan ukladu rownan
	  * @param	Vector &b	- wektor wyrazow wolnych ukladu rownan
	  *
	  */
	static void GaussJordanAlgorithmWithParalelization(Matrix& A, Vector& x, Vector& b);

	/**
  * Metoda rozszerzajaca metode GaussJordanAlgotithm() o dyrektywe
  * wektoryzujaca obliczenia
  *
  * @see GaussJordanAlgotithm()
  * @param	Matrix &A	- macierz wspolczynnikow ukladu rownan
  * @param	Vector &x	- wektor rozwiazan ukladu rownan
  * @param	Vector &b	- wektor wyrazow wolnych ukladu rownan
  *
  */
	static void GaussJordanAlgorithmWithVectorization(Matrix& A, Vector& x, Vector& b);
};

#endif // GAUSS_H