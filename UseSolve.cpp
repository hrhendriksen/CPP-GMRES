/*
* UseQ9.cpp
*
*  Created on: 25 apr. 2018
*      Author: hiddehendriksen
*      Use Q9 is linked to Q9.cpp via the Makefile:
*
*   all: UseQ9
*	Q9.o : Q9.cpp Q9.cpp
*		g++ -c -O Q9.cpp
*
*	UseQ9 : Q9.o UseQ9.cpp
*		g++ -O -o UseQ9 Q9.o UseQ9.cpp
*
*	might add print matrix function
*	remember: void function with pointer is dangerous, rather let function output pointer
*
 */

#include <iostream>
#include "Q9.hpp"

int main()
{
	double A[3][3] = { {1.0, 3.0, -1.0}, {1.0, 1.0, -1.0}, {3.0,11.0,15.0}}; // matrix A
	double B[3] = {9,1,35};		// vector b
	int n =3;
	double** a;
//	a = new double* [n];
//	for (int i = 0; i<n; i++)
//	{
//		a[i] = new double[n];
//	}
	a = allocatematrix(n,n);
	double* b = new double[3];

	for(int i=0; i<n; i++)
	{
		b[i] = B[i];

		for(int j=0; j<n; j++)
		{
			a[i][j] = A[i][j];
		}
	}
//	std::cout<<b[0];
	double** aug = create_aug(a, b, n);

	Gaussian_elimination(aug, n);

	solve_triangular(aug, n);
}



