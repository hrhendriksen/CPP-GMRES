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

	Vector c22(c_values,16);

	double d[4] = {2,4,5,2};
	Vector d2(d,4);

	print(a_matrix);
	std::cout<<"------------------------------ \n";
	// print();
	// std::cout<<"------------------------------ \n";
	// Matrix d = c_matrix + a_matrix;
	Matrix c2 = (c_matrix*c_matrix)/12;
	std::cout<<"---------------c=--------------- \n";
	print(c_matrix);
	std::cout<<"has dimensions"<< size(c_matrix)<<"\n";
	std::cout<<"---------------c22=--------------- \n";
	Vector ans = c_matrix*d2;
	std::cout <<"ans is "<< ans <<" with dimensions"<<size(ans)<<"\n";
	std::cout<<"---------------c33=--------------- \n";
	Vector ans2 = d2 * c_matrix;
	std::cout <<"ans is "<< ans2 <<"\n";
	std::cout<<"-------- I ---------------------- \n";
	print(diag(c22,0));
	print(diag(c22,1));
	print(diag(c22,3));

	std::cout<<"------------------------------ \n";
	// Matrix f(3,2);
	// try{
	// f = a_matrix-c_matrix;
	// }
	// catch(Exception& e){
	// 	e.DebugPrint();
	// }
	Vector b(d, 4);
	// double cc = b.Read(2);
	// std::cout<<cc<<"\n";
	Vector zz = c_matrix/b;
	// std::cout << gmres(c_matrix, b, zz, 10,10)<< "\n";
	return 0;
}
