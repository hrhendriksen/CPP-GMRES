#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "gmres.hpp"
#include <time.h>       /* time */

int main(int argc, char const *argv[])
{
	/* Testcase 1; example from the NLA course*/
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
	// Vector res = b_1-A_1*sol_1;
	// std::cout << "Test : "<<  norm(res) <<"\n";

	/* Testcase 2 - a less trivial 5x5 example*/
	// double A_2_values[25] =
	// {7,	3,	8,  2,	4,	 
	// 12,	13,	9,	10,	5,
	// 6,	11,	14, 15, 16,
	// 17,	18,	19,	20,	21,
	// 22,	23,	24,	25,	26};

	// Matrix A_2(A_2_values,5,5);

	// double x0arr_2[5] = {0,0,0,0,0};
	// Vector x0_2(x0arr_2,5);

	// double barr_2[5] = {1,0,3,-0.1,2};
	// Vector b_2(barr_2,5);

	// Vector sol_2 = gmres(A_2, b_2, x0_2, 100, 1e-5);	
	// std::cout << "The solution of the problem " << sol_2 << "\n";
	// Vector res = b_2-A_2*sol_2;
	// std::cout << "Test : "<<  norm(res) <<"\n";

	/* Testcase 3 - a random non sparse test example*/
	// int matrix_size;
	// std::cout<<"What matrix size do you want?\n";
	// std::cin >> matrix_size;

	// int max_iter;
	// std::cout<<"How many GMRES iterations, would you maximally like?\n";
	// std::cin >> max_iter;
	
	// double A_3_values[int(pow(matrix_size,2))];
	// double barr_3[matrix_size];
	// double x0arr_3[matrix_size];

	// srand(1);
	// for (int i = 0; i < int(pow(matrix_size,2)); ++i)
	// {
	// 	A_3_values[i] = rand()%10+1;
	// }	

	// for (int i = 0; i < matrix_size; ++i)
	// {
	// 	barr_3[i] = rand()%10+1;
	// 	x0arr_3[i] = rand()%10+1;
	// }

	// Matrix A_3(A_3_values, matrix_size, matrix_size);
	// // std::cout<<"========GMRES with A is ========\n";
	// // print(A_3);
	// Vector x0_3(x0arr_3, matrix_size);
	// // std::cout<<"========== x0 ============= is \n";
	// // std::cout<< x0_3 << "\n";
	// Vector b_3(barr_3, matrix_size);
	// // std::cout<<"========== b ============= is \n";
	// // std::cout<< b_3 << "\n";
	
	// Vector sol_3 = gmres(A_3, b_3, x0_3, max_iter+1, 1e-6);
	// std::cout << "The solution of the problem " << sol_3 << "\n";
	// Vector res = b_3-A_3*sol_3;
	// std::cout << "Test : "<<  norm(res) <<"\n";

	/* Testcase 4 - a sparse test case*/
	for (int test = 1; test < 201; ++test)
	{	int testsize =10*test;

		int max_iter = .2*testsize;
		double a[testsize-1];
		double b[testsize];
		double c[testsize-1];
		for (int i = 0; i < testsize; ++i)
		{
			a[std::max(0,i-1)] = rand()%10+1;
			b[i] =  rand()%10+1;
			c[std::max(0,i-1)] = rand()%10+1;
		}
		sparse_trid A_sp(testsize, a,b,c);

		Vector x0_sp(b,testsize);
		Vector b_sp(b,testsize);

		clock_t t1;
		t1 = clock();
		Vector sol_sp = gmres(A_sp, b_sp, x0_sp, max_iter+1, 1e-6);
		t1 = clock() - t1;
		std::cout<<"\t Sparse \t" << ((float)t1)/CLOCKS_PER_SEC<< "\t" << " n is " << testsize <<"\t";
		
		// std::cout << "The solution of the problem " << sol_sp << "\n";
		Vector res_sp = b_sp-A_sp*sol_sp;
		std::cout << "Sparse_Test : "<<  norm(res_sp) <<"\t";
//---------------------------------------------------------------
		/* Testcase 5 - a sparse test case*/
		// Copy the sparse matrix of Testcase 4 in to a Dense matrix.
		Matrix D(testsize, testsize);

		D(1,1) = b[0];
		D(1,2) = a[0];

		for (int i = 2; i < testsize; ++i)
		{
			for (int j = i-1; j <= i+1; ++j)
			{
				if (i-j == 0)
				{
					D(i,j) = b[i-1];
				}

				if (i-j == -1)
				{
					D(i,j) = a[i-1];
				}

				if (i-j == 1)
				{
					D(i,j) = c[i-1];
				}
			}
		}

		D(testsize,testsize-1) = c[testsize-1];
		D(testsize,testsize)   = b[testsize-1];
//---------------------------------------------------------------
		clock_t t2;
		t2 = clock();
		Vector sol_D = gmres(D, b_sp, x0_sp, max_iter+1, 1e-6);
		t2 = clock() - t2;
		std::cout<< "\t Dense \t" << ((float)t2)/CLOCKS_PER_SEC<< "\t";
		// std::cout << "The solution of the problem " << sol_sp << "\n";
		Vector res_D = b_sp-D*sol_D;
		std::cout << "Dense_Test : "<<  norm(res_sp) <<"\n";
	}

	return 0;
}