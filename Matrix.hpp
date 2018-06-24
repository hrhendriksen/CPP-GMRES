#ifndef MATRIXDEF
#define MATRIXDEF
//  Class written in such a way that code similar to Matlab code may be written
#include <cmath>
#include "Vector.hpp"
#include "Exception.hpp"
//  **********************
//  *  Class of matrices *
//  **********************

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
	Matrix(const Matrix& m); //copy constructor from a matrix
	Matrix(const Vector& v); //Copy (conversion) constructor from a vector
	
	//construct matrix of given size
	Matrix(int sizeVal_r,int sizeVal_c);
	
	// Overloaded initialisation constructor
	Matrix(double array[], int sizeVal_r, int sizeVal_c);
	Matrix(Vector& vector, int sizeVal_r, int sizeVal_c);
	
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
	friend Vector operator*(const Vector &v,  const Matrix& m);
	friend Vector operator*(const Matrix& m, const Vector &v);
	friend Matrix operator*(const Matrix& m, const double& a);
	friend Matrix operator/(const Matrix& m, const double& a);
	friend Matrix operator/(const Matrix& m, const double& a);

	//unary operator
	friend Matrix operator-(const Matrix& m);

	//other operators
	//Create augmented matrix
	friend Matrix create_aug(const Vector& v, const Matrix& m);
	// Reshape matrix
	friend Matrix cut(const Matrix& m, int new_m, int new_n);
	// Apply Gaussian Elimination
	void Gaussian_elimination();
	//Solve triangular system
	Vector solve_triangular();

	//Overloaded assignment
	Matrix& operator=(const Matrix& m);
	//indexing starting from 1
	double& operator()(int row, int column);
	//Friend functions
	friend void print(const Matrix& m);
	friend double det(const Matrix& m);
	friend std::pair<int,int> size(const Matrix& m);
	friend Vector operator/(const Matrix& m, const Vector& v);
	friend void writetoCSV(const Matrix& m);
	friend std::ostream& operator<<(std::ostream& output, const Matrix& m);
};

// All "friend" external operators and functions are declared as friend inside the class
// but their actual prototype definitions occur outside the class (here).
Matrix operator+(const Matrix& m1, const Matrix& m2);
Matrix operator-(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m1, const Matrix& m2);
Matrix operator*(const double& a,  const Matrix& m2);
Vector operator*(const Matrix& m,  const Vector &v);
Vector operator*(const Vector &v,  const Matrix& m);
Matrix operator*(const Matrix& m1, const double& a);
Matrix operator/(const Matrix& m1, const double& a);

//Create augmented matrix
Matrix create_aug(const Vector& v, const Matrix& m);
// Reshape matrix
Matrix cut(const Matrix& m, int new_m, int new_n);
// Matrix Gaussian_elimination(Matrix aug);
Vector operator/(const Matrix& m, const Vector& v);
//Unary operator
Matrix operator-(const Matrix& m);
// Declaration of size friend function
std::pair<int,int> size(const Matrix& m);
//Matrix output
std::ostream& operator<<(std::ostream& output, const Matrix& m);
// Write matrix friend function
void writetoCSV(const Matrix& m);
// Print matrix friend function
void print(const Matrix& m);
// Determinant friend function
double det(const Matrix& m);
// Function that generates n dimensional identity
Matrix eye(int n);
// Function that generates m by n zero matrix
Matrix zeros(int m, int n);
// Funciton that creates diagonal matrix with v on diagonal.
Matrix diag(const Vector& v,int n);

//  ********************************************
//  *  Class of sparse tridiagonal matrices *
//  ********************************************
class sparse_trid
{
private:
	//member variables
	int mSize;    // number of rows of the matrix
	double* superdiagonal;
	double* diagonal;
	double* subdiagonal;
public:
	// constructors
	// No default constructor
	//constructor sparse trid matrix of given size
	sparse_trid(int sizeVal);
	//constructor sparse trid matrix of with diagonals given as arrays
	sparse_trid(int sizeVal, double superd[], double d[], double subd[]);
	// Destructor
	~sparse_trid();
	//Function to get number of rows
	int GetNumberofRows() const;
	//Function to get number of columnss
	int GetNumberofColumns() const;
	// Print function
	friend void print(const sparse_trid& S);
	// Multiplication operator for a sparse tridiag matrix and a vector
	friend Vector operator*(const sparse_trid& S, Vector& v);
	// Function that converts a sparse matrix to a dense matrix
	friend Matrix sparse_trid2dense(const sparse_trid& S);
};

// Print function
void print(const sparse_trid& S);
// Multiplication friend function
Vector operator*(const sparse_trid& S, Vector& v);
// Conversion function from sparse to dense
Matrix sparse_trid2dense(const sparse_trid& S);
#endif