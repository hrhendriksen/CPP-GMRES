#ifndef VECTORDEF
#define VECTORDEF

//  **********************
//  *  Class of vectors  *
//  **********************

//  Class written in such a way that code similar to Matlab
//  code may be written

#include <cmath>
#include <utility>
#include "Exception.hpp"//  This class throws errors using the class "error"

class Vector
{
private:
   // member variables
   double* mData;   // data stored in vector
   int mSize;      // size of vector

public:
   // constructors
   // No default constructor
   // overridden copy constructor
   Vector(const Vector& v1);
   // construct vector of given size
   Vector(int sizeVal);
   // initialisation constructor
   Vector(double arr[], int sizeVal);

   // destructor
   ~Vector();

   // All "friend" external operators and functions are declared as friend inside the class (here)
   // but their actual prototype definitions occur outside the class.
   // Binary operators
   friend Vector operator+(const Vector& v1, const Vector& v2);
   friend Vector operator-(const Vector& v1, const Vector& v2);
   friend double operator*(const Vector& v1, const Vector& v2);
   friend Vector operator*(const Vector& v, const double& a);
   friend Vector operator*(const double& a, const Vector& v);
   friend Vector operator/(const Vector& v, const double& a);

   // Unary operator
   friend Vector operator-(const Vector& v);

   //other operators
   //assignment
   Vector& operator=(const Vector& v);

   //indexing
   double& operator()(int i);

   //reading
   const double& operator()(int i) const;
   
   //norm (as a member method)
   double norm(int p=2) const;

   // Read vector entry
   double Read(int i) const;

   // functions that are friends
   //output
   friend std::ostream& operator<<(std::ostream& output, const Vector& v);
   // Get vector length
   friend int length(const Vector& v);
   // Get vector norm
   friend double norm(const Vector& v, int p);
   // Get vector size
   friend std::pair<int,int> size(const Vector& v);
   // Cut off the vector
   friend Vector cut(const Vector& v, int new_m);
   // Write the vector into a CSV output file
   friend void writetoCSV(const Vector& v);
};

// All "friend" external operators and functions are declared as friend inside the class
// but their actual prototype definitions occur outside the class (here).
// Binary operators
Vector operator+(const Vector& v1, const Vector& v2);
Vector operator-(const Vector& v1, const Vector& v2);
double operator*(const Vector& v1, const Vector& v2);
Vector operator*(const Vector& v, const double& a);
Vector operator*(const double& a, const Vector& v);
Vector operator/(const Vector& v, const double& a);

// Unary operator
Vector operator-(const Vector& v);

// function prototypes
// Prototype signature of length() friend function
int length(const Vector& v);
// Prototype signature of norm() friend function
double norm(const Vector& v, int p=2);
// Prototype signature of size() friend function
std::pair<int,int> size(const Vector& v);
// Prototype signature of cut() friend function
Vector cut(const Vector& v, int new_m);
// Prototype signature of writetoCSV() friend function
void writetoCSV(const Vector& v);

#endif
