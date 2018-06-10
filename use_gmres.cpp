#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "gmres.hpp"

int main(int argc, char const *argv[])
{
	double c_values[25] = {0, 0, 0, 0, 1,
   						   1, 0, 0, 0, 1,
   						   0, 1, 0, 0, 0,
   						   0, 0, 1, 0, 0, 
   						   0, 0, 0, 1, 0};

	Matrix c_matrix(c_values,5,5);
	
	double x0a[5] = {0,0,0,0,0};
	Vector x0(x0a,5);

	double ba[5] = {1,0,0,0,0};
	Vector b(ba,5);
	std::cout << gmres(c_matrix, b, x0, 6, .1)<< "\n";
	return 0;
}