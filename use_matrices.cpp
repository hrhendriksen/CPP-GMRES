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
   	double c_values[16] = {2, 4, 3, 1,
   						   1, 1, 0, 2,
   						   0, 1, 1, 0,
   						   3, 0, 1, 2};

	Matrix c_matrix(c_values,4,4);


	a_matrix.print();
	std::cout<<"------------------------------ \n";
	c_matrix.print();
	std::cout<<"------------------------------ \n";
	// Matrix d = c_matrix + a_matrix;
	Matrix c2 = (c_matrix*c_matrix)/12;
	std::cout<<"---------------c*c=--------------- \n";
	c2.print();
	std::cout<<"---------------c*c=--------------- \n";
	// d.print();
	// Matrix e = -c_matrix;
	// std::cout<<"------------------------------ \n";
	// // e.print();
	// std::cout<<"------------------------------ \n";
	// Matrix f(3,2);
	// try{
	// f = a_matrix-c_matrix;
	// }
	// catch(Exception& e){
	// 	e.DebugPrint();
	// }
	double d[4] = {2,4,5,2};
	Vector b(d, 4);
	// double cc = b.Read(2);
	// std::cout<<cc<<"\n";
	Vector zz = c_matrix/b;
	std::cout << zz<< "\n";
	return 0;
}
