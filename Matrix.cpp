#include <iostream>
#include "Matrix.hpp"
#include <fstream>
//construct sparse trid matrix of given size
// the first element of the superdiagonal is set to zero but is not part of the matrix,
// the last element of the subdiagonal is set to zero but is not part of the matrix.
sparse_trid::sparse_trid(int sizeVal)
{
	mSize = sizeVal;

	superdiagonal = new double [sizeVal];
	diagonal = new double [sizeVal];
	subdiagonal = new double [sizeVal];

	for (int i = 0; i < sizeVal; ++i)
	{
		superdiagonal[i] = 0.0;
		diagonal[i] = 0.0;
		subdiagonal[i] = 0.0;
	}
}

// Construct sparse trid matrix with 3 arrays
sparse_trid::sparse_trid(int sizeVal, double superd[], double d[], double subd[])
{
	mSize = sizeVal;

	superdiagonal = new double [sizeVal];
	diagonal = new double [sizeVal];
	subdiagonal = new double [sizeVal];

	superdiagonal[0] = 0.0;
	diagonal[0] = d[0];
	subdiagonal[0] = subd[0];

	for (int i = 1; i < sizeVal-1; ++i)
	{
		superdiagonal[i] = superd[i-1];
		diagonal[i]=d[i];
		subdiagonal[i] = subd[i];
	}
	superdiagonal[sizeVal-1] = superd[sizeVal-2];
	diagonal[sizeVal-1] = d[sizeVal-1];
	subdiagonal[sizeVal-1] = 0.0;
}

sparse_trid::~sparse_trid()
{
	delete[] superdiagonal;
	delete[] diagonal;
	delete[] subdiagonal;
}

// Function to get number of rows of matrix
int sparse_trid::GetNumberofRows() const
{
	return mSize;
}

void print(const sparse_trid& S)
{
	std::cout<<"superdiagonal: \n";
	std::cout<<"[";
	for (int i = 0; i < S.mSize; ++i)
	{
		std::cout<<S.superdiagonal[i];
		if (i<S.mSize-1)
		{
			std::cout<< ", ";
		}
	}
	std::cout<<"]\n";

	std::cout<<"diagonal: \n";
	std::cout<<"[";
	for (int i = 0; i < S.mSize; ++i)
	{
		std::cout<<S.diagonal[i];
		if (i<S.mSize-1)
		{
			std::cout<< ", ";
		}
	}
	std::cout<<"]\n";

	std::cout<<"subdiagonal: \n";
	std::cout<<"[";
	for (int i = 0; i < S.mSize; ++i)
	{
		std::cout<<S.subdiagonal[i];
		if (i<S.mSize-1)
		{
			std::cout<< ", ";
		}
	}
	std::cout<<"]\n";
}

Vector operator*(const sparse_trid& S, Vector& v)
{
	assert(S.mSize == length(v));
	Vector product(S.mSize);
	product(1) = S.diagonal[0]*v(1)+S.superdiagonal[1]*v(2);
	for (int i = 2; i < S.mSize; ++i)
	{	
		product(i)=v(i-1)*S.subdiagonal[i-2]+v(i)*S.diagonal[i-1]+v(i+1)*S.superdiagonal[i];
	}
	product(S.mSize) = S.diagonal[S.mSize-1]*v(S.mSize)+S.subdiagonal[S.mSize-2]*v(S.mSize-1);
	return product;
}

//copy constructor
Matrix::Matrix(const Matrix& m)
{
	mSize_r = m.mSize_r;
	mSize_c = m.mSize_c;

	mData = new double* [m.mSize_r];

	for(int i = 0; i<m.mSize_r; i++)
	{
		mData[i] = new double [mSize_c];
		for (int j = 0; j < mSize_c; ++j)
		{
			mData[i][j] = m.mData[i][j];
		}
	}
}
// copy (conversion) constructor that makes a Matrix out of a Vector
Matrix::Matrix(const Vector& v)
{	
	int m = length(v);
	mSize_r = m;
	mSize_c = 1;

	mData = new double* [m];
	for (int i = 0; i < m; ++i)
	{
		mData[i] = new double [1];
		
		mData[i][0] = v.Read(i+1);
	}
}

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

// initialisation constructor
Matrix::Matrix(Vector& vector, int sizeVal_r, int sizeVal_c)
{
		mSize_r = sizeVal_r;
		mSize_c = sizeVal_c;
		mData = new double* [sizeVal_r];

		for(int i = 0; i<sizeVal_r; i++)
		{
			mData[i] = new double [sizeVal_c];
			for (int j = 0; j < sizeVal_c; ++j)
			{
				mData[i][j] = vector.Read(sizeVal_c*i+j+1);
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

// Function to get number of rows of matrix
int Matrix::GetNumberofRows() const
{
	return mSize_r;
}

//Function to get number of columns of matrix
int Matrix::GetNumberofColumns() const
{
	return mSize_c;
}


// definition of + between two matrices
Matrix operator+(const Matrix& m1, const Matrix& m2)
{
	int m;
	int n;

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

std::ostream& operator<<(std::ostream& output, const Matrix& m) 
{	
  for (int i=0; i< m.mSize_r; i++)
    { 
    	for (int j = 0; j < m.mSize_c; ++j)
    	{
    		output << m.mData[i][j];
    		if (j != m.mSize_c-1)
			{
				output  << "\t";
    		}
      		else
			{
				output  << std::endl;
			}
    	}
    }
  return output;  // for multiple << operators.
}

void writetoCSV(const Matrix& m)
{
	std::ofstream out("output.csv");
	assert(out.is_open());
	out << m;
	out.close();
	return;
}



// Function that multiplies two matrices
Matrix operator*(const Matrix& m1, const Matrix& m2)
{
	Matrix B(m1.mSize_r,m2.mSize_c);
	assert(m1.mSize_c == m2.mSize_r);

	for (int i = 0; i < m1.mSize_r; ++i)
	{
		for (int j = 0; j < m2.mSize_c; ++j)
		{
			B.mData[i][j]=0.0;
			for (int k = 0; k < m1.mSize_c; ++k)
			{
				B.mData[i][j] += B.mData[i][j] + m1.mData[i][k]*m2.mData[k][j];
			}
		}
	}
	return B;
}

// Function that multiplies matrix and scalar
Matrix operator*(const double& a, const Matrix& m2)
{
	for (int i = 0; i < m2.mSize_r; ++i)
	{
		for (int j = 0; j < m2.mSize_c; ++j)
		{
			m2.mData[i][j] = a*m2.mData[i][j];
		}
	}
	return m2;
}

// Function that multiplies scalar and matrix
Matrix operator*(const Matrix& m2, const double& a)
{
	return a*m2;
}

Matrix operator/(const Matrix& m1, const double& a)
{
  if (a == 0.0)
  {
     throw Exception("div 0", "Attempt to divide by zero");
  }

  return m1*(1.0/a); 
}

Vector operator*(const Matrix& m,  const Vector& v)
{
	assert(m.mSize_c == length(v));
	double array[m.mSize_r];

	for (int i = 0; i < m.mSize_r; ++i)
	{	
		array[i] = 0;
		for (int j = 0; j < m.mSize_c; ++j)
		{
			array[i] += m.mData[i][j]*v.Read(j+1);	
		}
	}
	Vector ans(array, m.mSize_r);
	return ans;
}

Vector operator*(const Vector &v, const Matrix& m)
{
	assert(m.mSize_r == length(v));
	double array[m.mSize_c];
	for (int j = 0; j < m.mSize_c; ++j)
	{
		array[j]=0;

		for (int i = 0; i < m.mSize_r; ++i)
		{
			array[j] += m.mData[i][j]*v.Read(i+1);
		}
	}
	Vector ans(array, m.mSize_c);
	return ans;
}

// Function that creates the augmented matrix from a vector and a matrix
Matrix create_aug(const Vector& v, const Matrix& m)
{
	assert(m.mSize_r == length(v));
	Matrix aug(m.mSize_r, m.mSize_c+1);

	for(int i = 1; i<m.mSize_r+1; i++)
	{
		for(int j = 1; j<m.mSize_c+2; j++)
		{
			if(j<m.mSize_c+1){
					aug(i,j) = m.mData[i-1][j-1];
					// std::cout<<aug(i,j)<<"\t";
					}
			else if(j == m.mSize_c+1)
				{
				aug(i,j) = v.Read(i);
				// std::cout<<aug(i,j)<<"\t";
				}
		}
			// std::cout<<"\n";
	}
	return aug;
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

// return the size of a matrix
Vector size(const Matrix& m)
{
	double values[2] = {m.mSize_r,m.mSize_c};
	return Vector(values,2);
}


// print matrix
void print(const Matrix& m1)
{
	for (int i = 0; i < m1.mSize_r; ++i)
	{

		for (int j = 0; j < m1.mSize_c; ++j)
		{
			std::cout<< m1.mData[i][j] << "\t";
		}
	std::cout<<"\n";
	}
}

// Recursive function to calculate determinant
double det(const Matrix& m)
{
	assert(m.mSize_r == m.mSize_c);
	double determinant = 0.0;
	
	if (m.mSize_r == 1)
	{
		determinant = m.mData[0][0];
	}
	else
	{
		// More than one entry of matrix
		for (int i_outer = 0; i_outer < m.mSize_r; ++i_outer)
			{
				Matrix sub_matrix(m.mSize_r-1, m.mSize_c-1);
				for (int i = 0; i < m.mSize_r-1; ++i)
				{
					for (int j = 0; j < i_outer; ++j)
					{
							sub_matrix(i+1,j+1) = m.mData[i+1][j];
					}
					for (int j = i_outer; j < m.mSize_c-1; ++j)
					{
						sub_matrix(i+1,j+1) = m.mData[i+1][j+1];
					}
				}
			double sub_matrix_determinant = det(sub_matrix);

			determinant += pow(-1.0,i_outer) * m.mData[0][i_outer]*sub_matrix_determinant;
			}	
	}
	return determinant;
}

// Cut matrix
Matrix cut(const Matrix& m, int new_m, int new_n)
{
	int min_rows;
	int min_columns;

	if (m.mSize_r < new_m)
	{
		min_rows = m.mSize_r;
	}
	else
	{
		min_rows = new_m;
	}

	if(m.mSize_c < new_n)
	{
		min_columns = m.mSize_c;
	}
	else
	{
		min_columns = new_n;
	}
	
	Matrix new_matrix(new_m,new_n);
	for (int i = 0; i < min_rows; ++i)
	{
		for (int j = 0; j < min_columns; ++j)
		{
			new_matrix.mData[i][j] = m.mData[i][j];
		}
	}

	return new_matrix;
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

Vector operator/(const Matrix& m, const Vector& v)
{
	assert(m.mSize_c == m.mSize_r);

	if (det(m)!= 0)
	{
		throw Exception("inverse non-existent", "Matrix is non-singular");	
	}

	Matrix augmented = create_aug(v, m);
	Matrix check = augmented.Gaussian_elimination();

	return check.solve_triangular();	
}

Matrix Matrix::Gaussian_elimination()
{	
	int m = mSize_r;
	int n = mSize_c;

	// std::cout << "Gaussian_elimination started with m is"<< m <<"\n";
	// for(int pp=0; pp<m; pp++)
	// 	{
	// 		for (int kk = 0; kk < n; ++kk)
	// 		{
	// 			std::cout << mData[pp][kk] << "\t";
	// 		}
	// 		std::cout << "\n";
	// 	}
	int h = 0;
	int k = 0;
	while( h<m and k<m)
	{
		// find pivot
		int max_p=h;
		for(int l = h; l<m; l++)
		{
			if(std::fabs(mData[l][k])>std::fabs(mData[max_p][k]))
			{
				max_p = l;
			}
		}
	//		std::cout << max_p;
			// Now we know what the max pivot is.
			// swap rows
	//		std::cout << "here";
		if(max_p != h)
		{
			double ph[m+1];
			for(int kk = 0; kk<m+1; kk++)
				{
					ph[kk] =mData[h][kk];
				}

			for(int jj = 0; jj<m+1; jj++)
			{
				mData[h][jj] = mData[max_p][jj];
			}

			for(int ll = 0; ll<m+1; ll++)
			{
				mData[max_p][ll] = ph[ll];
			}
		}
			// std::cout << "___________________ \n";
	
			// for(int pp=0; pp<m; pp++)
			// {
			// 	for (int kk = 0; kk < n; ++kk)
			// 	{
			// 		std::cout << mData[pp][kk] << "\t";
			// 	}
			// 	std::cout << "\n";
			// }


		if(mData[h][k]!=0)
		{
			for(int mm = (h+1); mm <m; mm++)
			{
				double AA = (mData[mm][k]/mData[h][k]);
				for(int oo = k; oo<m+1; oo++)
				{
					mData[mm][oo] = mData[mm][oo] - AA*mData[h][oo];
				}
			}
		}
		k+=1;
		h+=1;
	}

	double *array = new double[m*n];
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			array[i*n+j]=mData[i][j];
		}
	}
	return Matrix(array, m, n);
}

Vector Matrix::solve_triangular()
{
	// std::cout<<"solve_triangular\n";
	int n = mSize_c;
	int m = mSize_r;
	// std::cout<<"number of columns: "<< n << "\n";
	// std::cout<<"number of rows: "<< m << "\n";

	// std::cout<<"(2,2) entry"<<GE.mData[2][2]<<"\n";
	// std::cout<<"(2,3) entry"<<GE.mData[2][3]<<"\n";

	double x[m];
	for (int i=m-1; i>=0; i--)
	{
		x[i] = mData[i][n-1]/mData[i][i];
		// std::cout<<"running value x[i]"<<x[i]<<"\n";

		for (int k=i-1;k>=0; k--)
		{
			mData[k][n-1] -= mData[k][i] * x[i];
		}
		// std::cout << GE.mData[1][3]
	}
	// std::cout<<"_._._._._._._. \n";
	for(int fin = 0; fin<m; fin++)
	{
		// std::cout << x[fin]<<"\n";
		// std::cout<<"print vector";
	 }
	return Vector(x,m);
}

Matrix eye(int n)
{
	double data[n*n];
	for (int i = 0; i < n*n; ++i)
	{
		if (i%(n+1) == 0)
		{
			data[i] = 1;
		}
		else
		{
			data[i]=0;
		}
	}
	Matrix I(data,n,n);
	return I;
}

Matrix zeros(int m, int n)
{
	Matrix M(m,n);
	return M;
}

Matrix diag(const Vector& v, int k)
{	int n = length(v)+std::abs(k);
	std::cout << "n is " << n << "\n";
	std::cout << "length of v is"<< length(v)<<"\n";
	double data[n*n];
	int l = 0;
	for (int i = 0; i < n*n; ++i)
	{
		if ((i-k)%(n+1) == 0 && l < length(v))
		{	
			std::cout << "i = " << i << "\t l = " << l << "\n";
			data[i] = v.Read(l+1);
			l+=1;
		}

		else
		{
			data[i]=0;
		}
	}
	Matrix D(data,n,n);
	return D;
}


