#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"
#include <chrono>

int main(int argc, char* argv[])
{	
	// Test constructors
	// Test initialisation from array
   	double a_values[16] = {2, 4, 3, 1,
   						   1, 1, 0, 2,
   						   0, 1, 1, 0,
   						   3, 0, 1, 2};
	Matrix a_matrix(a_values,4,4);
	print(a_matrix);
	// Test initialisation from size integers
	Matrix b_matrix(2,3);
	print(b_matrix);

	// Test copy constructor from matrix;
	Matrix c_matrix(a_matrix);
	print(c_matrix);
	 // Test copy constructor from vector
	Vector b(a_values, 16);
	Matrix kk(b);
	print(kk);

	// Test exception
  	std::cout << "The following throws an exception\n";
	Matrix BB(4,4);
	try
	{
	BB = c_matrix/0.0;
	}

  	catch (Exception &ex)
	{
      ex.DebugPrint();
  	}
	double d[4] = {2,4,5,2};
	Vector d2(d,4);
	std::cout<<"-+----------------------------- \n";
	Matrix d_matrix = c_matrix + a_matrix - 2*a_matrix;
	print(d_matrix);
	std::cout<<"-*/----------------------------- \n";
	Matrix c2 = (c_matrix*c_matrix)/12;
	print(c2);
	
	//Test size
	std::cout<<"-size----------------------------- \n";
	std::cout<<"c has dimensions ("<< size(c_matrix).first << ", " << size(c_matrix).second <<")\n";
	
	//Test MVP
	std::cout<<"---------------MVP--------------- \n";
	Vector ans = c_matrix*d2;
	std::cout <<"ans is "<< ans <<"\n";
	std::cout<<"---------------VMP=--------------- \n";
	Vector ans2 = d2 * c_matrix;
	std::cout <<"ans is "<< ans2 <<"\n";

	// Test Identity
	std::cout<<"-------- I ---------------------- \n";
	print(eye(2));
	
	// Test diag
	std::cout<<"-------- diag ---------------------- \n";
	print(diag(d2,0));
	print(diag(d2,1));
	print(diag(d2,3));
	std::cout<<"------------------------------ \n";
	Matrix f(4,4);
	try{
	f = a_matrix-c_matrix;
	}
	catch(Exception& e){
		e.DebugPrint();
	}
	
	// Test cut function
	std::cout<<"---------------- cut -------------- \n";
	Matrix c_cut_matrix = cut(c_matrix, 5, 5);
	print(c_cut_matrix);

	//Test determinant and backslash
	std::cout<<"-------------- det and backslash ---------- \n";

	std::cout<< "det(c) is :"<< det(c_matrix) << "\n";
	Vector zz = c_matrix/d2;
	std::cout << "the solution is:"<< zz << "\n";

	//Test sparse_trid class
	std::cout<<"-------------- sparse trid and sparse trid mvp ---------- \n";
	sparse_trid KKK(3);
	double az[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	double bz[17] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
	double cz[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	sparse_trid LLL(17,az,bz,cz);
	print(LLL);

	Vector LVL(bz,17);
	Vector product = LLL*LVL;
	std::cout << "MVP" << product << "\n";

/* ======================================================
 				Testcase 1 - Simple MVP.
======================================================*/
// Create random sparse tridiagonal matrix, ts is testsize
/* initialize random seed: */
// srand (time(NULL));
// for (int ts = 2; ts <= 3000; ++ts)
// 	{
// 		// Create sparse matrix
// 		double super[ts-1], diag[ts], sub[ts-1], ran[ts];
// 		for (int i = 0; i < ts; ++i)
// 		{
// 			super[std::max(0,i-1)] = rand()%10+1;
// 			diag[i] =  rand()%10+1;
// 			sub[std::max(0,i-1)] = rand()%10+1;
// 			ran[i] = rand()%10+1;
// 		}
// 		sparse_trid A_sp(ts, super,diag,sub);
// 		// Create random vector
// 		Vector RAN(ran, ts);
// 		// Copy into a dense matrix
// 		Matrix D = sparse_trid2dense(A_sp);
// 		// Clock both mvps and output times in table
// 		auto start1 = std::chrono::high_resolution_clock::now();
// 		Vector sparse = A_sp*RAN;
// 		auto end1 = std::chrono::high_resolution_clock::now();
// 		auto diff1 = end1 - start1;
// 		std::cout << " n is  " << ts <<"\t";
// 		std::cout<<"\t Sparse \t" << std::chrono::duration 
// 		<double, std::milli> (diff1).count() << " \t ms" << "\t";
// 		auto start2 = std::chrono::high_resolution_clock::now();
// 		Vector dense = D*RAN;
// 		auto end2 = std::chrono::high_resolution_clock::now();
// 		auto diff2 = end2 - start2;
// 		std::cout<<"\t Dense \t" << std::chrono::duration <double, std::milli> (diff2).count() << "\t ms" << "\n";
// 	}
	return 0;
}
