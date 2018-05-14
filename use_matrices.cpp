#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"

int main(int argc, char* argv[])
{	
	double a_values[2] = {2,1};
	Matrix a_matrix(a_values,2,1);
	// Matrix b_matrix(2,3);
   	double c_values[9] = {1000, 2, 3, 17, 50, 10, 1, 2, 3};
	Matrix c_matrix(c_values,3,3);


	a_matrix.print();
	std::cout<<"------------------------------ \n";
	c_matrix.print();
	std::cout<<"------------------------------ \n";
	// Matrix d = c_matrix + a_matrix;
	// d.print();
	// Matrix e = -c_matrix;
	std::cout<<"------------------------------ \n";
	// e.print();
	std::cout<<"------------------------------ \n";
	Matrix f(3,2);
	try{
	f = a_matrix-c_matrix;
	}
	catch(Exception& e){
		e.DebugPrint();
	}
	double d[3] = {1,2,3};
	Vector b(d, 3);
	create_aug(b, c_matrix).print();

	return 0;
}