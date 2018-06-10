#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "gmres.hpp"

int main(int argc, char const *argv[])
{
	/* Testcase 1*/
	// double A_1_values[25] = 
	// {0, 0, 0, 0, 1,
 //   	 1, 0, 0, 0, 0,
 //   	 0, 1, 0, 0, 0,
 //   	 0, 0, 1, 0, 0, 
 //   	 0, 0, 0, 1, 0};

	// Matrix A_1(A_1_values,5,5);
	
	// double x0arr_1[5] = {0,0,0,0,0};
	// Vector x0_1(x0arr_1,5);

	// double barr_1[5] = {1,0,0,0,0};
	// Vector b_1(barr_1,5);
	// Vector sol_1 = gmres(A_1, b_1, x0_1, 6, .1);
	// std::cout << "The solution of the problem " << sol_1 << "\n";
	// std::cout << "Test : "<< b_1-A_1*sol_1 <<"\n";

	/* Testcase 2*/
	double A_2_values[25] =
	{7,	3,	8,  2,	4,	 
	12,	13,	9,	10,	5,
	6,	11,	14, 15, 16,
	17,	18,	19,	20,	21,
	22,	23,	24,	25,	26};

	Matrix A_2(A_2_values,5,5);

	double x0arr_2[5] = {0,0,0,0,0};
	Vector x0_2(x0arr_2,5);

	double barr_2[5] = {1,0,3,-0.1,2};
	Vector b_2(barr_2,5);

	Vector sol_2 = gmres(A_2, b_2, x0_2, 10, .1);	
	std::cout << "The solution of the problem " << sol_2 << "\n";
	std::cout << "Test : "<< b_2-A_2*sol_2 <<"\n";

	return 0;
}