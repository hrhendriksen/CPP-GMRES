#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "gmres.hpp"

int main(int argc, char const *argv[])
{
	double c_values[16] = {2, 4, 3, 1,
   						   1, 1, 0, 2,
   						   0, 1, 1, 0,
   						   3, 0, 1, 2};

	Matrix c_matrix(c_values,4,4);
	double d[4] = {2,4,5,2};
	Vector d2(d,4);
	std::cout << gmres(c_matrix, d2, d2, 10, .1)<< "\n";
	return 0;
}