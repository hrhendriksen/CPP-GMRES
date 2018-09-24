#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "gmres.hpp"
#include <chrono>
#include <time.h>       /* time */

int main(int argc, char const *argv[])
{
	/* ====================================================
 				Testcase 1 - example from the NLA course.
	======================================================*/
 
	double A_1_values[25] = 
	{0, 0, 0, 0, 1,
   	 1, 0, 0, 0, 0,
   	 0, 1, 0, 0, 0,
   	 0, 0, 1, 0, 0, 
   	 0, 0, 0, 1, 0};

	Matrix A_1(A_1_values,5,5);
	
	double x0arr_1[5] = {0,0,0,0,0};
	Vector x0_1(x0arr_1,5);

	double barr_1[5] = {1,0,0,0,0};
	Vector b_1(barr_1,5);
	Vector sol_1 = gmres(A_1, b_1, x0_1, 6, .1);
	Vector res = b_1-A_1*sol_1;
	std::cout << "Test : "<<  norm(res) <<"\n";

/* ====================================================
 				Testcase 2 - a less trivial 5x5 example.
	======================================================*/
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

	Vector sol_2 = gmres(A_2, b_2, x0_2, 100, 1e-5);	
	std::cout << "The solution of the problem " << sol_2 << "\n";
	Vector res_2 = b_2-A_2*sol_2;
	std::cout << "Test : "<<  norm(res_2) <<"\n";

/* ====================================================
 			Testcase 3 - a random non sparse test example.
	======================================================*/
	int matrix_size = 600;
	int max_iter = 1*matrix_size+1;
	int max_adds  = 100;

	// Preallocate matrix to store residuals in 
	Matrix Res_mat(max_iter, max_adds);
	for (int j = 1; j <= max_adds; ++j)
	{

		double A_3_values[int(pow(matrix_size,2))];
		double barr_3[matrix_size];
		double x0arr_3[matrix_size];

		// Create random matrix and vector
		srand(1);
		for (int ii = 0; ii < int(pow(matrix_size,2)); ++ii)
		{
			A_3_values[ii] = rand()%10+1;
		}	

		for (int iii = 0; iii < matrix_size; ++iii)
		{
			barr_3[iii] = rand()%10+1;
			x0arr_3[iii] = rand()%10+1;
		}

		Matrix A_3(A_3_values, matrix_size, matrix_size);
		std::cout<<"========GMRES with A is ========\n";
		print(A_3);
		Vector x0_3(x0arr_3, matrix_size);
		std::cout<<"========== x0 ============= is \n";
		std::cout<< x0_3 << "\n";
		Vector b_3(barr_3, matrix_size);
		std::cout<<"========== b ============= is \n";
		std::cout<< b_3 << "\n";
		
		//Let GMRES output the vector of residuals
		Matrix ZZZ = A_3+j*eye(matrix_size);
		Vector sol_3 = gmres(ZZZ, b_3, x0_3, max_iter+1, 1e-6);

		//Write residual into matrix
		for (int i = 1; i < length(sol_3); ++i)
		{
			Res_mat(i,j) = sol_3(i);
		}
		
	}

	writetoCSV(Res_mat);


/* ====================================================
 			Testcase 3 - Diagonally dominant example.
	======================================================*/
		int Ttestsize = 10000;
		int Mmax_iter = .6*Ttestsize;

		double Aa[Ttestsize-1];
		double Bb[Ttestsize];
		double Cc[Ttestsize-1];
		double Xx[Ttestsize];

		for (int i = 0; i < Ttestsize; ++i)
		{
			Aa[std::max(0,i-1)] = -1;
			Bb[i] =  3;
			Cc[std::max(0,i-1)] = -1.1;
			Xx[i]=rand()%2;
		}

		sparse_trid Aa_sp(Ttestsize, Aa,Bb,Cc);
		// print(Aa_sp);
		Vector Xxx(Xx, Ttestsize);
		Vector Bbb = Aa_sp*Xxx;

		time_t Tt;
		Tt = clock();
		Vector Ssol = gmres(Aa_sp, Bbb, Bbb, Mmax_iter+1, 1e-6);
		Tt = clock() - Tt;
		std::cout<< "\t Sparse \t" << ((float)Tt)/CLOCKS_PER_SEC<< "\t";
		// std::cout << "The solution of the problem " << Ssol << "\n";
		Vector Rres = Bbb-Aa_sp*Ssol;
		std::cout << "Sparse_Test : "<<  norm(Rres) <<"\n";

	/* Testcase 5a - a sparse test case*/
	// initialize random seed:
	srand (time(NULL));
	int number_of_tests = 100;
	for (int test = 1; test <= number_of_tests; ++test)
	{	
		int testsize =10*test;
		int max_iter = testsize+1;
		double a[testsize-1], b[testsize], c[testsize-1];
//====NDD==============================================================		
		// for (int i = 0; i < testsize; ++i)
		// {
		// 	a[std::max(0,i-1)] = rand()%10+1;
		// 	b[i] =  rand()%10+1;
		// 	c[std::max(0,i-1)] = rand()%10+1;
		// }
//=====DD==============================================================
		for (int i = 0; i < testsize; ++i)
		{
			a[std::max(0,i-1)] = rand()%10+1;
			b[i] = 50 +(rand()%10+1);
			c[std::max(0,i-1)] = rand()%10+1;
		}
//=====================================================================
		sparse_trid A_sp(testsize, a,b,c);
		Vector x0_sp(b,testsize);
		Vector b_sp(b,testsize);
		auto start1 = std::chrono::high_resolution_clock::now();
		Vector sol_sp = gmres(A_sp, b_sp, x0_sp, max_iter+1, 1e-6);
		auto end1 = std::chrono::high_resolution_clock::now();
		auto diff1 = end1 - start1;
		std::cout<<"Sparse \t" << std::chrono::duration 
		<double, std::milli> (diff1).count() << " \t ms" << "\t";
		Vector res_sp = b_sp-A_sp*sol_sp;
		std::cout << "S_Test : "<<  norm(res_sp)/norm(b_sp) <<"\t";
//---------------------------------------------------------------
		Matrix D = sparse_trid2dense(A_sp);
		auto start2 = std::chrono::high_resolution_clock::now();
		Vector sol_D = gmres(D, b_sp, x0_sp, max_iter+1, 1e-6);
		auto end2 = std::chrono::high_resolution_clock::now();
		auto diff2 = end2 - start2;
		std::cout<<"\t Dense \t" << std::chrono::duration <double, std::milli> (diff2).count() << "\t ms";
		Vector res_D = (b_sp-D*sol_D);
		std::cout << "D_T : "<<  norm(res_sp)/norm(b_sp) <<"\n";
	}

	return 0;
}