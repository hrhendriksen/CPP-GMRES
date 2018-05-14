/*
 * Q9.cpp
 *
 *  Created on: 24 apr. 2018
 *      Author: hiddehendriksen
 */
#include <iostream>
#include "Q9.hpp"
#include <cmath>
#include "Matrix.hpp"
// #include "vector.hpp"
Matrix create_aug(Matrix& A, Vector& b)
{
	int n = A.mSize_c;
	// double** aug = allocatematrix(n, n+1);
	// Matrix aug(n,n);
	// for(int i = 0; i<n; i++){
	// 	for(int j = 0; j<n+1; j++){
	// 		if(j<n){
	// 			aug[i][j] = A[i][j];
	// 			std::cout<<aug[i][j]<<"\t";
	// 		}
	// 		else if(j == n)
	// 		{
	// 			aug[i][j] = b[i];
	// 			std::cout<<aug[i][j]<<"\t";
	// 		}
	// 	}
	// 	std::cout<<"\n";
	// }
	// return aug;
	std::cout<<n;
}

// double** Gaussian_elimination(double** aug, int n)
// {
// 	std::cout<<"GE started \n";
// 	for(int pp=0; pp<n; pp++)
// 		{
// 			std::cout << aug[pp][0] << "\t" << aug[pp][1]<<"\t"<< aug[pp][2]<< "\t"<< aug[pp][3]<<"\n";
// 		}
// 	int h = 0;
// 	int k = 0;
// 	while( h<n and k<n){
// 		// find pivot
// 		int max_p=h;
// 		for(int l = h; l<n; l++)
// 		{
// 			if(std::fabs(aug[l][k])>std::fabs(aug[max_p][k]))
// 			{
// 				max_p = l;
// 			}
// 			}
// 	//		std::cout << max_p;
// 			// Now we know what the max pivot is.
// 			// swap rows
// 	//		std::cout << "here";
// 			if(max_p != h)
// 			{
// 				double ph[n+1];
// 				for(int kk = 0; kk<n+1; kk++)
// 					{
// 						ph[kk] =aug[h][kk];
// 					}
// //			std::cout << ph[0] << "\t" << ph[1]<<"\t"<< ph[2]<< "\t" << ph[3]<<"\n";
// 				for(int jj = 0; jj<n+1; jj++)
// 				{
// 					aug[h][jj] = aug[max_p][jj];
// 				}

// 				for(int ll = 0; ll<n+1; ll++)
// 				{
// 					aug[max_p][ll] = ph[ll];
// 				}
// 			}
// 			std::cout << "___________________ \n";
// 			for(int pp=0; pp<n; pp++)
// 			{
// 				std::cout << aug[pp][0] << "\t" << aug[pp][1]<<"\t"<< aug[pp][2]<< "\t"<< aug[pp][3]<<"\n";
// 			}


// 			if(aug[h][k]!=0)
// 			{
// 				for(int mm = (h+1); mm <n; mm++)
// 				{
// 					double AA = (aug[mm][k]/aug[h][k]);
// 					for(int oo = k; oo<n+1; oo++)
// 					{
// 	//					std::cout << "o = " << oo;
// 						aug[mm][oo] = aug[mm][oo] - AA*aug[h][oo];
// 					}
// 				}
// 			}
// 			std::cout << "___________________ \n";
// 			for(int pp=0; pp<n; pp++)
// 			{
// 				std::cout << aug[pp][0] << "\t" << aug[pp][1]<<"\t"<< aug[pp][2]<< "\t"<< aug[pp][3]<<"\n";
// 			}

// 			k+=1;
// 			h+=1;
// 		}
// 	return aug;
// }
// void solve_triangular(double** aug, int n)
// {
// 	double x[n];
// 		for (int i=2; i>=0; i--)
// 		{
// 			x[i] = aug[i][n]/aug[i][i];

// 			for (int k=i-1;k>=0; k--)
// 			{
// 				aug[k][n] -= aug[k][i] * x[i];
// 			}

// 		}
// 		std::cout<<"_._._._._._._. \n";
// 		for(int fin = 0; fin<n; fin++)
// 		{
// 			std::cout << x[fin]<<"\n";
// 		}
// }

// double** allocatematrix(int size_a, int size_b)
// {
// 	double**MA;
// 	MA = new double* [size_a];
// 	for (int i = 0; i<size_a; i++)
// 	{
// 		MA[i] = new double[size_b];
// 	}
// 	return MA;
// }



