#ifndef MATRIXDEF
#define MATRIXDEF

//  **********************
//  *  Class of vectors  *
//  **********************

//  Class written in such a way that code similar to Matlab code may be written
#include <cmath>
#include "Vector.hpp"
#include "Exception.hpp"

class Matrix
{
private: 
	// member variables
	double** mData; // data stored in a matrix (entries)
	int mSize_r;    // number of rows of the matrix
	int mSize_c;	// number of columns of the matrix
public:
	//constructors
	// No default constructor
	
	// overriden copy constructor
	Matrix(const Matrix& m1); //copy constructor
	
	//construct matrix of given size
	Matrix(int sizeVal_r,int sizeVal_c);
	
	// initialisation constructor
	Matrix(double arr[], int sizeVal_r, int sizeVal_c);
	
	//destructor
	~Matrix();	

	//Function to get number of rows
	int GetNumberofRows() const;

	//Function to get number of columns	
	int GetNumberofColumns() const;


   // All "friend" external operators and functions are declared as friend inside the class (here)
   // but their actual prototype definitions occur outside the class.
   // Binary operators

	friend Matrix operator+(const Matrix& m1, const Matrix& m2);
	friend Matrix operator-(const Matrix& m1, const Matrix& m2);
	friend Matrix operator*(const Matrix& m1, const Matrix& m2);
	friend Matrix operator*(const double& a,  const Matrix& m2);
	friend Matrix operator*(const Matrix& m1, const double& a);
	friend Matrix operator/(const Matrix& m1, const double& a);
	// friend Matrix operator/(const Matrix& m1, const Vector& v1);

	//unary operator
	friend Matrix operator-(const Matrix& m);
	
	//other operators
	//Create augmented matrix
	friend Matrix create_aug(const Vector& v, const Matrix& m);
	// Apply Gaussian Elimination
	Matrix Gaussian_elimination();
	//Solve triangular system
	Vector solve_triangular();


	//Overloaded assignment
	Matrix& operator=(const Matrix& m);
	//indexing starting from 1
	double& operator()(int row, int column);
	//output
	void print();
	// friend std::ostream& operator<<(std::ostream& output, const Matrix& m);

};

// All "friend" external operators and functions are declared as friend inside the class
// but their actual prototype definitions occur outside the class (here).
// Binary operators
Matrix operator+(const Matrix& m1, const Matrix& m2);
Matrix operator-(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m1, const Matrix& m2);
Matrix operator*(const double& a,  const Matrix& m2);
Matrix operator*(const Matrix& m1, const double& a);
Matrix operator/(const Matrix& m1, const double& a);
Matrix create_aug(const Vector& v, const Matrix& m);
// Matrix Gaussian_elimination(Matrix aug);
Vector operator/(const Matrix& m, const Vector& v);
// //Unary operator
Matrix operator-(const Matrix& m);

// void print(const Matrix& m);


#endif