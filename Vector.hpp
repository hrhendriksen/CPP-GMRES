#ifndef VECTORDEF
#define VECTORDEF

//  **********************
//  *  Class of vectors  *
//  **********************


//  Class written in such a way that code similar to Matlab
//  code may be written

#include <cmath>
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

   // Get vector length

   // Read vector entry
   double Read(int i) const;



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
   friend int length(const Vector& v);


   //other operators
   //assignment
   Vector& operator=(const Vector& v);
   //indexing
   double& operator()(int i);
   //output
   friend std::ostream& operator<<(std::ostream& output, const Vector& v);

   //norm (as a member method)
   double norm(int p=2) const;
   // functions that are friends
   friend double norm(Vector& v, int p);
   friend Vector size(const Vector& v);
   friend Vector cut(const Vector& v, int new_m);
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
Vector cut(const Vector& v, int new_m);
int length(const Vector& v);

// Unary operator
Vector operator-(const Vector& v);

// function prototypes
double norm(Vector& v, int p=2);
// Prototype signature of length() friend function
int length(const Vector& v);
Vector size(const Vector& v);

void writetoCSV(const Vector& v);
#endif
