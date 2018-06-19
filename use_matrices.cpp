#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"
#include <chrono>

int main(int argc, char* argv[])
{	
	// double a_values[2] = {2,1};
	// Matrix a_matrix(a_values,2,1);
	// // Matrix b_matrix(2,3);
 //   	double c_values[16] = {2, 4, 3, 1,
 //   						   1, 1, 0, 2,
 //   						   0, 1, 1, 0,
 //   						   3, 0, 1, 2};

	// Matrix c_matrix(c_values,4,4);
	// double zero = 0;
	// Matrix BB(4,4);
	// BB = c_matrix/zero;
	// Vector c22(c_values,16);

	// double d[4] = {2,4,5,2};
	// Vector d2(d,4);

	// print(c_matrix);
	// std::cout<<"------------------------------ \n";
	// print();
	// std::cout<<"------------------------------ \n";
	// Matrix d = c_matrix + a_matrix;
	// Matrix c2 = (c_matrix*c_matrix)/12;
	// std::cout<<"---------------c=--------------- \n";
	// print(c_matrix);
	// std::cout<<"has dimensions"<< size(c_matrix)<<"\n";
	// std::cout<<"---------------c22=--------------- \n";
	// Vector ans = c_matrix*d2;
	// std::cout <<"ans is "<< ans <<" with dimensions"<<size(ans)<<"\n";
	// std::cout<<"---------------c33=--------------- \n";
	// Vector ans2 = d2 * c_matrix;
	// std::cout <<"ans is "<< ans2 <<"\n";
	// std::cout<<"-------- I ---------------------- \n";
	// print(diag(c22,0));
	// print(diag(c22,1));
	// print(diag(c22,3));
	// std::cout<<"------------------------------ \n";
	// Matrix f(3,2);
	// try{
	// f = a_matrix-c_matrix;
	// }
	// catch(Exception& e){
	// 	e.DebugPrint();
	// }
	// Vector b(d, 4);
	// Matrix kk(b);
	// // std::cout<<"----------------KK-------------- \n";
	// // print(cut(c_matrix,5,5));
	// // std::cout<<"------------------------------ \n";
	// // double cc = b.Read(2);
	// // std::cout<<cc<<"\n";
	// // Vector zz = c_matrix/b;
	// // std::cout << gmres(c_matrix, d2, d2, 10,10)<< "\n";
	// sparse_trid KKK(3);
	// double az[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	// double bz[17] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
	// double cz[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	// sparse_trid LLL(17,az,bz,cz);
	// // print(LLL);
	// Vector LVL(bz,17);
	// Vector product = LLL*LVL;
	// // std::cout << product << "\n";

/* Testcase 1 - Simple MVP.*/
// Create random sparse tridiagonal matrix, ts is testsize
for (int ts = 2; ts <= 3000; ++ts)
	{
		/* initialize random seed: */
  		srand (time(NULL));
		double super[ts-1], diag[ts], sub[ts-1], ran[ts];
		for (int i = 0; i < ts; ++i)
		{
			super[std::max(0,i-1)] = rand()%10+1;
			diag[i] =  rand()%10+1;
			sub[std::max(0,i-1)] = rand()%10+1;
			ran[i] = rand()%10+1;
		}
		sparse_trid A_sp(ts, super,diag,sub);
		Vector RAN(ran, ts);
		// Copy into a dense matrix
		Matrix D = sparse_trid2dense(A_sp);
		auto start1 = std::chrono::high_resolution_clock::now();
		Vector sparse = A_sp*RAN;
		auto end1 = std::chrono::high_resolution_clock::now();
		auto diff1 = end1 - start1;
		std::cout << " n is  " << ts <<"\t";
		std::cout<<"\t Sparse \t" << std::chrono::duration 
		<double, std::milli> (diff1).count() << " \t ms" << "\t";
		auto start2 = std::chrono::high_resolution_clock::now();
		Vector dense = D*RAN;
		auto end2 = std::chrono::high_resolution_clock::now();
		auto diff2 = end2 - start2;
		std::cout<<"\t Dense \t" << std::chrono::duration <double, std::milli> (diff2).count() << "\t ms" << "\n";
	}
	return 0;
}
// 