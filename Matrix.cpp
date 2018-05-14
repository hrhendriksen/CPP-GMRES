#include <iostream>
#include "Matrix.hpp"

// constructor that creates vector of given size with
// double precision entries all initially set to zero

Matrix::Matrix(int sizeVal_r, int sizeVal_c)
{
	mSize_r = sizeVal_r;
	mSize_c = sizeVal_c;
	// std::cout<<mSize_r;

	mData = new double* [sizeVal_r];

	for(int i = 0; i<sizeVal_r; i++)
	{
		mData[i] = new double [sizeVal_c];
		for (int j = 0; j < sizeVal_c; ++j)
		{
			mData[i][j] = 0.0;
		}
	}
}

//copy constructor
Matrix::Matrix(const Matrix& m1)
{
	mSize_r = m1.mSize_r;
	mSize_c = m1.mSize_c;

	mData = new double* [m1.mSize_r];

	for(int i = 0; i<m1.mSize_r; i++)
	{
		mData[i] = new double [mSize_c];
		for (int j = 0; j < mSize_c; ++j)
		{
			mData[i][j] = m1.mData[i][j];
		}
	}
}

// initialisation constructor
Matrix::Matrix(double arr[], int sizeVal_r, int sizeVal_c)
{
	mSize_r = sizeVal_r;
	mSize_c = sizeVal_c;

	mData = new double* [sizeVal_r];

	for(int i = 0; i<sizeVal_r; i++)
	{
		mData[i] = new double [sizeVal_c];
		for (int j = 0; j < sizeVal_c; ++j)
		{
			mData[i][j] = arr[sizeVal_c*i+j];
		}
	}
}

// Matrix desctructor
Matrix::~Matrix()
{
	for (int i = 0; i < mSize_r; ++i)
	{
		delete[] mData[i];
	}
	delete[] mData;
}
// print matrix
void Matrix::print()
{
	for (int i = 0; i < mSize_r; ++i)
	{

		for (int j = 0; j < mSize_c; ++j)
		{
			std::cout<< mData[i][j] << "\t";
		}
	std::cout<<"\n";
	}
}

// definition of + between two matrices
Matrix operator+(const Matrix& m1, const Matrix& m2)
{
	int m;
	int n;
	// std::cout<< m1.mSize_r <<m1.mSize_c;
	// std::cout<< m2.mSize_r <<m2.mSize_c;
	if (m1.mSize_r > m2.mSize_r)
	{
		m = m1.mSize_r; 
	}
	else
	{
		m = m2.mSize_r;
	}

	if (m1.mSize_c > m2.mSize_c)
	{
		n = m1.mSize_c; 
	}
	else
	{
		n = m2.mSize_c;
	}

	Matrix w(m,n);
	std::cout<<m<<n;
	// add the matrices
	if (m1.mSize_r == m2.mSize_r) 
	{
		// std::cout<<m1.mSize_c<<m2.mSize_c;
		if (m1.mSize_c == m2.mSize_c)
		{
			std::cout<<"equaltty";
			for (int i = 0; i < m1.mSize_r; ++i)
			{
				for(int j=0; j<m1.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
					}
			}	
	 	}
		else if(m1.mSize_c < m2.mSize_c)
		{

			for (int i = 0; i < m1.mSize_r; ++i)
				{
					for (int j = 0; j < m1.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
					}

					for (int j=m1.mSize_c; j < m2.mSize_c; j++)
					{
						w.mData[i][j] = m2.mData[i][j];
					}
				}
				std::cerr<<"matrix add - different number of colums\n";
				std::cerr<<"extra entries of smaller matrix assumed to be 0.0\n";
		}
		else //m1.mSize_c > m2.mSize_c
		{
			for (int i = 0; i < m1.mSize_r; ++i)
			{
				for (int j = 0; j < m2.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
					}	
				for (int j = m2.mSize_c; j < m1.mSize_c; ++j)
				{
					w.mData[i][j] = m1.mData[i][j];
				}
			}
				std::cerr<<"matrix add (WHY)- different number of colums\n";
				std::cerr<<"extra entries of smaller matrix assumed to be 0.0\n";
		}
	}
	//different number of rows:	
	else if (m1.mSize_r < m2.mSize_r)
	{
		if (m1.mSize_c == m2.mSize_c)
		{
			for (int j = 0; j < m1.mSize_c; ++j)
			{
				for(int i=0; i<m1.mSize_r; ++i)
					{
					w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
					}
				for (int i = m1.mSize_r; i < m2.mSize_r; ++i)
				{
					w.mData[i][j] = m2.mData[i][j];
				}
			}	
		}
		else if(m1.mSize_c < m2.mSize_c)
		{

			for (int i = 0; i < m1.mSize_r; ++i)
				{
					for (int j = 0; j < m1.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
					}
				}
			for (int i = m1.mSize_r; i < m2.mSize_r; ++i)
			{					
				for (int j = 0; j < m2.mSize_c; j++)
				{
					w.mData[i][j] = m2.mData[i][j];
				}
			}

			for (int j = m1.mSize_c; j < m2.mSize_c; ++j)
			{
				for (int i = 0; i < m2.mSize_r; ++i)
				{
					w.mData[i][j] = m2.mData[i][j];
				}
			}
		
			std::cerr<<"matrix add - different number of colums\n";
			std::cerr<<"extra entries of smaller matrix assumed to be 0.0\n";

		
		}
		else //m1.mSize_c > m2.mSize_c and m1.mSize_r < m2.mSize_r
		{
			for (int i = 0; i < m1.mSize_r; ++i)
			{
				for (int j = 0; j < m2.mSize_c; ++j)
				{
					w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
				}

				for (int j = 0; j < m1.mSize_r; ++j)
				{
					w.mData[i][j] = m1.mData[i][i];
				}
			}
			for (int i = m1.mSize_r; i < m2.mSize_r; ++i)
			{
				for (int j = 0; i < m2.mSize_c; ++j)
				{
					w.mData[i][j] = m2.mData[i][j];
				}
			}
		}
		std::cerr<<"bbbmatrix add - different number of colums\n";
		std::cerr<<"extra entries of smaller matrix assumed to be 0.0\n";
	}
	else // m1.mSize_r > m2.mSize_r
	{
	if (m1.mSize_c == m2.mSize_c)
		{
			for (int j = 0; j < m2.mSize_c; ++j)
			{
				for(int i=0; i<m2.mSize_r; ++i)
					{
					w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
					}
				for (int i = m2.mSize_r; i < m1.mSize_r; ++i)
				{
					w.mData[i][j] = m1.mData[i][j];
				}
			}	
		}
		else if(m2.mSize_c < m1.mSize_c)
		{

			for (int i = 0; i < m2.mSize_r; ++i)
				{
					for (int j = 0; j < m2.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
					}
				}
			for (int i = m2.mSize_r; i < m1.mSize_r; ++i)
			{					
				for (int j=0; j < m1.mSize_c; j++)
				{
					w.mData[i][j] = m1.mData[i][j];
				}
			}

			for (int j = m2.mSize_c; j<m1.mSize_c; ++j)
			{
				for (int i = 0; i < m1.mSize_r; ++i)
				{
					w.mData[i][j] = m1.mData[i][j];
				}
			}
		
			std::cerr<< "matrix add - different number of colums\n";
			std::cerr<< "extra entries of smaller matrix assumed to be 0.0\n";

		
		}
		else //m1.mSize_c < m2.mSize_c and m1.mSize_r > m2.mSize_r
		{
			for (int i = 0; i < m2.mSize_r; ++i)
			{
				for (int j = 0; j < m1.mSize_c; ++j)
				{
					w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
				}
			}

			for (int j = 0; j < m2.mSize_c; ++j)
			{
				for (int i = m2.mSize_r; i < m1.mSize_r; ++i)
				{
					w.mData[i][j] = m1.mData[i][j];	
				}
				
			}
			
			for (int i = 0; i < m2.mSize_r; ++i)
			{
				for (int j = m1.mSize_c; j < m2.mSize_c; ++j)
				{
					w.mData[i][j] = m2.mData[i][j];
				}
			}

		std::cerr<< "???matrix add - different number of colums\n";
		std::cerr<< "extra entries of smaller matrix assumed to be 0.0\n";
		}

	}
	return w;
}

//definition of the unary operator - 
Matrix operator-(const Matrix& m)
{
	Matrix w(m.mSize_r, m.mSize_c);

	for (int i = 0; i < m.mSize_r; ++i)
	{
		for (int j = 0; j < m.mSize_c; ++j)
		{
			w.mData[i][j] = -m.mData[i][j];
		}
	}
	return w;
}

//definition of the - minus operator
Matrix operator-(const Matrix& m1, const Matrix& m2)
{
	int m;
	int n;
	// std::cout<< m1.mSize_r <<m1.mSize_c;
	// std::cout<< m2.mSize_r <<m2.mSize_c;
	if (m1.mSize_r > m2.mSize_r)
	{
		m = m1.mSize_r; 
	}
	else
	{
		m = m2.mSize_r;
	}

	if (m1.mSize_c > m2.mSize_c)
	{
		n = m1.mSize_c; 
	}
	else
	{
		n = m2.mSize_c;
	}

	Matrix w(m,n);
	std::cout<<m<<n;
	// add the matrices
	if (m1.mSize_r == m2.mSize_r) 
	{
		// std::cout<<m1.mSize_c<<m2.mSize_c;
		if (m1.mSize_c == m2.mSize_c)
		{
			std::cout<<"equaltty";
			for (int i = 0; i < m1.mSize_r; ++i)
			{
				for(int j=0; j<m1.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]+m2.mData[i][j];
					}
			}	
	 	}
		else if(m1.mSize_c < m2.mSize_c)
		{

			for (int i = 0; i < m1.mSize_r; ++i)
				{
					for (int j = 0; j < m1.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
					}

					for (int j=m1.mSize_c; j < m2.mSize_c; j++)
					{
						w.mData[i][j] = -m2.mData[i][j];
					}
				}
				std::cerr<<"matrix add - different number of colums\n";
				std::cerr<<"extra entries of smaller matrix assumed to be 0.0\n";
		}
		else //m1.mSize_c > m2.mSize_c
		{
			for (int i = 0; i < m1.mSize_r; ++i)
			{
				for (int j = 0; j < m2.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
					}	
				for (int j = m2.mSize_c; j < m1.mSize_c; ++j)
				{
					w.mData[i][j] = m1.mData[i][j];
				}
			}
				std::cerr<<"matrix add (WHY)- different number of colums\n";
				std::cerr<<"extra entries of smaller matrix assumed to be 0.0\n";
		}
	}
	//different number of rows:	
	else if (m1.mSize_r < m2.mSize_r)
	{
		if (m1.mSize_c == m2.mSize_c)
		{
			for (int j = 0; j < m1.mSize_c; ++j)
			{
				for(int i=0; i<m1.mSize_r; ++i)
					{
					w.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
					}
				for (int i = m1.mSize_r; i < m2.mSize_r; ++i)
				{
					w.mData[i][j] = -m2.mData[i][j];
				}
			}	
		}
		else if(m1.mSize_c < m2.mSize_c)
		{

			for (int i = 0; i < m1.mSize_r; ++i)
				{
					for (int j = 0; j < m1.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
					}
				}
			for (int i = m1.mSize_r; i < m2.mSize_r; ++i)
			{					
				for (int j = 0; j < m2.mSize_c; j++)
				{
					w.mData[i][j] = -m2.mData[i][j];
				}
			}

			for (int j = m1.mSize_c; j < m2.mSize_c; ++j)
			{
				for (int i = 0; i < m2.mSize_r; ++i)
				{
					w.mData[i][j] = -m2.mData[i][j];
				}
			}
		
			std::cerr<<"matrix add - different number of colums\n";
			std::cerr<<"extra entries of smaller matrix assumed to be 0.0\n";

		
		}
		else //m1.mSize_c > m2.mSize_c and m1.mSize_r < m2.mSize_r
		{
			for (int i = 0; i < m1.mSize_r; ++i)
			{
				for (int j = 0; j < m2.mSize_c; ++j)
				{
					w.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
				}

				for (int j = 0; j < m1.mSize_r; ++j)
				{
					w.mData[i][j] = m1.mData[i][i];
				}
			}
			for (int i = m1.mSize_r; i < m2.mSize_r; ++i)
			{
				for (int j = 0; i < m2.mSize_c; ++j)
				{
					w.mData[i][j] = -m2.mData[i][j];
				}
			}
		}
		std::cerr<<"bbbmatrix add - different number of colums\n";
		std::cerr<<"extra entries of smaller matrix assumed to be 0.0\n";
	}
	else // m1.mSize_r > m2.mSize_r
	{
	if (m1.mSize_c == m2.mSize_c)
		{
			for (int j = 0; j < m2.mSize_c; ++j)
			{
				for(int i=0; i<m2.mSize_r; ++i)
					{
					w.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
					}
				for (int i = m2.mSize_r; i < m1.mSize_r; ++i)
				{
					w.mData[i][j] = m1.mData[i][j];
				}
			}	
		}
		else if(m2.mSize_c < m1.mSize_c)
		{

			for (int i = 0; i < m2.mSize_r; ++i)
				{
					for (int j = 0; j < m2.mSize_c; ++j)
					{
						w.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
					}
				}
			for (int i = m2.mSize_r; i < m1.mSize_r; ++i)
			{					
				for (int j=0; j < m1.mSize_c; j++)
				{
					w.mData[i][j] = m1.mData[i][j];
				}
			}

			for (int j = m2.mSize_c; j<m1.mSize_c; ++j)
			{
				for (int i = 0; i < m1.mSize_r; ++i)
				{
					w.mData[i][j] = m1.mData[i][j];
				}
			}
		
			std::cerr<< "matrix add - different number of colums\n";
			std::cerr<< "extra entries of smaller matrix assumed to be 0.0\n";

		
		}
		else //m1.mSize_c < m2.mSize_c and m1.mSize_r > m2.mSize_r
		{
			for (int i = 0; i < m2.mSize_r; ++i)
			{
				for (int j = 0; j < m1.mSize_c; ++j)
				{
					w.mData[i][j] = m1.mData[i][j]-m2.mData[i][j];
				}
			}

			for (int j = 0; j < m2.mSize_c; ++j)
			{
				for (int i = m2.mSize_r; i < m1.mSize_r; ++i)
				{
					w.mData[i][j] = m1.mData[i][j];	
				}
				
			}
			
			for (int i = 0; i < m2.mSize_r; ++i)
			{
				for (int j = m1.mSize_c; j < m2.mSize_c; ++j)
				{
					w.mData[i][j] = -m2.mData[i][j];
				}
			}

		std::cerr<< "???matrix add - different number of colums\n";
		std::cerr<< "extra entries of smaller matrix assumed to be 0.0\n";
		}

	}
	return w;
}

//definition of the matrix operation =
Matrix& Matrix::operator=(const Matrix& m)
{
// Check if matrices have the same dimensions,
// if rhs has smaller dimensions than lhs, assume missing entries are 0
// if rhs has larger dimensions than lhs, then throw
	if ((m.mSize_c > mSize_c)||(m.mSize_r> mSize_r))
	{
		throw Exception("dimensions mismatch", "matrix assignment operator - matrices have different dimensions");
	}
	else if ((m.mSize_r < mSize_r) && (m.mSize_c < mSize_c))
	{

		// Copy the rhs values in the lhs matrix
		for (int i = 0; i < m.mSize_r; ++i)
		{
			for (int j = 0; j < m.mSize_c; ++j)
			{
				mData[i][j] = m.mData[i][j];
			}
		}
		// Fill the remainder of the lhs matrix with zeroes
		for (int i = m.mSize_r; i < mSize_r; ++i)
		{
			for (int j = 0; j < mSize_c; ++j)
			{
				mData[i][j] = 0.0;	
			}
		}

		for (int j = m.mSize_c; j < mSize_c; ++j)
		{
			for (int i = 0; i < m.mSize_r; ++i)
			{
				mData[i][j] = 0.0;
			}
		}

	}
	return *this;
}

double& Matrix::operator()(int i, int j)
{
	if ((i < 1)||(j<1))
	{
		throw Exception("out of range", "accessing matrix through () - index too Small");
	}
	else if ((i>mSize_r)||(j>mSize_c))
	{
		throw Exception("out of range","accessing vector through () - index too high");
	}
	return mData[i-1][j-1];
}

Vector operator/(const Vector& v, const Matrix& m)
{
	Matrix aug = create_aug(v,m);
	return v;
}


Matrix create_aug(const Vector& v, const Matrix& m)
{
	assert(m.mSize_c == m.mSize_r);
	Matrix aug(m.mSize_r, m.mSize_c+1);

	// for(int i = 0; i<n; i++)
	// {
	// 	for(int j = 0; j<n+1; j++)
	// 	{
	// 		if(j<n){
	// 				aug(i,j) = m.mData[i][j];
	// 				std::cout<<aug(i,j)<<"\t";
	// 				}
	// 		else if(j == n)
	// 			{
	// 			aug[i][j] = v.mData[i];
	// 			std::cout<<aug[i][j]<<"\t";
	// 			}
	// 	}
	// 		std::cout<<"\n";
	// }
	return aug;
}


