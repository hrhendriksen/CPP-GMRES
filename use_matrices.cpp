#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Matrix.hpp"
#include <time.h>

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
	Vector b(d, 4);
	Matrix kk(b);
	// std::cout<<"----------------KK-------------- \n";
	// print(cut(c_matrix,5,5));
	// std::cout<<"------------------------------ \n";
	// double cc = b.Read(2);
	// std::cout<<cc<<"\n";
	// Vector zz = c_matrix/b;
	// std::cout << gmres(c_matrix, d2, d2, 10,10)<< "\n";
	sparse_trid KKK(3);
	double az[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	double bz[17] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
	double cz[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	sparse_trid LLL(17,az,bz,cz);
	// print(LLL);
	Vector LVL(bz,17);
	Vector product = LLL*LVL;
	// std::cout << product << "\n";

	// Testcase 1 - Simple Matrix Test.

	// Create random sparse tridiagonal matrix


for (int i = 2; i <= 3000; ++i)
	{
		int testsize = i;
		double aaa[testsize-1];
		double bbb[testsize];
		double ccc[testsize-1];
		double ddd[testsize];

		for (int i = 0; i < testsize; ++i)
		{
			aaa[std::max(0,i-1)] = rand()%10+1;
			bbb[i] =  rand()%10+1;
			ccc[std::max(0,i-1)] = rand()%10+1;
			ddd[i] = rand()%10+1;
		}

		sparse_trid A_sp(testsize, aaa,bbb,ccc);

		// Copy into a dense matrix
		Matrix D(testsize, testsize);
		Vector dDd(ddd, testsize);

		D(1,1) = bbb[0];
		D(1,2) = aaa[0];

		for (int i = 2; i < testsize; ++i)
		{
			for (int j = i-1; j <= i+1; ++j)
			{
				if (i-j == 0)
				{
					D(i,j) = bbb[i-1];
				}

				if (i-j == -1)
				{
					D(i,j) = aaa[i-1];
				}

				if (i-j == 1)
				{
					D(i,j) = ccc[i-1];
				}
			}
		}

		D(testsize,testsize-1) = ccc[testsize-1];
		D(testsize,testsize)   = bbb[testsize-1];

		clock_t t1;
		t1 = clock();
		Vector sparse = A_sp*dDd;
		t1 = clock() - t1;
		std::cout<<"\t Sparse \t" << ((float)t1)/CLOCKS_PER_SEC<< "\t" << " n is " << testsize <<"\t";


		clock_t t2;
		t2 = clock();
		Vector dense =  D*dDd;
		t2 = clock() - t2;
		std::cout<< "\t Dense \t" << ((float)t2)/CLOCKS_PER_SEC<< "\n";
	}
		
	return 0;
}
