#include <iostream>
#include <fstream>
#include "Vector.hpp"
#include "utility"

// copy constructor - creates vector with the same entries as v1
Vector::Vector(const Vector& v1)
{
  mSize = v1.mSize;
  mData=new double[mSize];
  for (int i=0; i<v1.mSize; i++)
  {
    mData[i] = v1.mData[i];
  }
}
// constructor that creates vector of given size with
// double precision entries all initially set to zero
Vector::Vector(int sizeVal)
{
  mData=new double[sizeVal];
  mSize = sizeVal;
  for (int i=0; i<mSize; i++)
  {
     mData[i] = 0.0;
  }
}

// initialisation constructor
Vector::Vector(double arr[], int sizeVal)
{
  mSize = sizeVal;
  mData=new double[mSize];
  for (int i=0; i<sizeVal; i++)
  {
    mData[i] = arr[i];
  }
}

// destructor - deletes pointer
Vector::~Vector()
{
  delete[] mData;
}



// definition of + between two vectors
Vector operator+(const Vector& v1, 
                        const Vector& v2)
{
  int n;

  //  set n to be the length of the longest vector and create a vector
  //  of that length to be returned
  if (v1.mSize > v2.mSize)
  {
      n = v1.mSize;
  }
  else
  {
      n = v2.mSize;
  }
  Vector w(n);


  //  add the vectors
  //  if one vector is shorter than the other assume missing entries are 0
  if (v1.mSize == v2.mSize)
  {
      for (int i=0; i<v1.mSize; i++)
       {
       w.mData[i] = v1.mData[i] + v2.mData[i];
  }
  }
  else if (v1.mSize > v2.mSize)
  {
      for (int i=0; i<v2.mSize; i++)
  {
    w.mData[i] = v1.mData[i] + v2.mData[i];
  }
      for (int i=v2.mSize; i<v1.mSize; i++)
  {
    w.mData[i] = v1.mData[i];
  }
      std::cerr<<"vector add - vectors different lengths\n";
      std::cerr<<"extra entries of shorter vector assumed to be 0.0\n";
  }

  else
  {
      for (int i=0; i<v1.mSize; i++)
  {
    w.mData[i] = v1.mData[i] + v2.mData[i];
  }
      for (int i=v1.mSize; i<v2.mSize; i++)
  {
    w.mData[i] = v2.mData[i];
  }
      std::cerr<<"vector add - vectors different lengths\n";
      std::cerr<<"extra entries of shorter vector assumed to be 0.0\n";
  }
  return w;
}

// definition of - between two vectors
Vector operator-(const Vector& v1, const Vector& v2)
{
  int n;

  //  set n to be the length of the longest vector and create a vector
  //  of that length to be returned
  if (v1.mSize > v2.mSize)
  {
    n = v1.mSize;
  }
  else
  {
    n = v2.mSize;
  }
  Vector w(n);

  //  subtract the vectors
  //  if one vector is shorter than the other assume missing entries are 0
  if (v1.mSize == v2.mSize)
  {
    for (int i=0; i<v1.mSize; i++)
    {
      w.mData[i] = v1.mData[i] - v2.mData[i];
    }
  }

  else if (v1.mSize > v2.mSize)
  {
    for (int i=0; i<v2.mSize; i++)
    {
      w.mData[i] = v1.mData[i] - v2.mData[i];
    }
    for (int i=v2.mSize; i<v1.mSize; i++)
    {
      w.mData[i] = v1.mData[i];
    }
    std::cerr<<"vector subtract - vectors different lengths\n";
    std::cerr<<"extra entries of shorter vector assumed to be 0.0\n";
  }

  else
  {
    for (int i=0; i<v1.mSize; i++)
    {
      w.mData[i] = v1.mData[i] - v2.mData[i];
    }
    for (int i=v1.mSize; i<v2.mSize; i++)
    {
      w.mData[i] = v2.mData[i];
    }
    std::cerr<<"vector subtract - vectors different lengths\n";
    std::cerr<<"extra entries of shorter vector assumed to be 0.0\n";
  }
  return w;
}

// definition of scalar product between two vectors
double operator*(const Vector& v1, const Vector& v2)
{
  int n;
  double dp;
  //  check vectors are of the same length
  //  if not assume missing entries are 0 
  if (v1.mSize > v2.mSize)
  {
    n = v2.mSize;
    std::cerr<<"scalar product - vectors different lengths\n";
    std::cerr<<"extra entries of shorter vector assumed to be 0.0\n";
  }
  else if (v1.mSize < v2.mSize)
  {
    n = v1.mSize;
    std::cerr<<"scalar product - vectors different lengths\n";
    std::cerr<<"extra entries of shorter vector assumed to be 0.0\n";
  }
  else
  {
    n = v1.mSize;
  }
  //  compute scalar product
  dp = 0.0;
  for (int i=0; i<n; i++)
  {
    dp += v1.mData[i] * v2.mData[i];
  }
  return dp;
}

// definition of multiplication between a vector and a scalar
Vector operator*(const Vector& v, const double& a)
{
//  create a vector of the same length as v with entries equal to a*v

  Vector w(v.mSize);

  for (int i=0; i<v.mSize; i++)
  {
    w.mData[i] = a * v.mData[i];
  }

  return w;
}

// definition of multiplication between a scalar and a vector
Vector operator*(const double& a, const Vector& v)
{
//  create a vector of the same length as v with entries equal to a*v
  Vector w(v.mSize);
  for (int i=0; i<v.mSize; i++)
  {
      w.mData[i] = a * v.mData[i];
  }
  return w;
}

// definition of division of a vector by a scalar
Vector operator/(const Vector& v, const double& a)
{
  if (a == 0.0)
  {
     throw Exception("div 0", "Attempt to divide by zero");
  }
//  create a vector of the same length as v with entries equal to v/a

  Vector w(v.mSize);
  for (int i=0; i<v.mSize; i++)
    {
      w.mData[i] = v.mData[i] / (double)a;
    }
  return w;
}


// definition of the unary operator -
Vector operator-(const Vector& v)
{
//  create a vector w with entries equal to -v
  Vector w(v.mSize);
  for (int i=0; i<v.mSize; i++)
  {
      w.mData[i] = -v.mData[i];
  }
  return w;
}

// definition of vector operator =
Vector& Vector::operator=(const Vector& v)
{
  //  check both vectors have same length
  //  if rhs vector is too short, assume missing entries are 0
  //  if rhs vector is too long then throw
  if (v.mSize > mSize)
  {
    throw Exception("length mismatch",
    "vector assignment operator - vectors have different lengths");
  }
  else if (v.mSize < mSize)
  {
    for (int i=0; i<v.mSize; i++)
    {
      mData[i] = v.mData[i];
    }
    for (int i=v.mSize; i<mSize; i++)
    {
      mData[i] = 0.0;
    }
    std::cerr << "vector assignment operator - copied vector was too short";
    std::cerr << " and has been extended with zeroes\n";
  }
  else
  {
    for (int i=0; i<mSize; i++)
    {
      mData[i] = v.mData[i]; 
    }
  }
  return *this;
}

// definition of vector operator ()
// allows v.mData[i] to be written as v(i+1), as in Matlab and FORTRAN
double& Vector::operator()(int i)
{
  if (i < 1)
  {
    throw Exception("out of range",
    "accessing vector through () - index too small");
  }
  else if (i > mSize)
  {
    throw Exception("length mismatch",
    "accessing vector through () - index too high");
  }
  return mData[i-1];
}

const double& Vector::operator()(int i) const
{
    if (i < 1)
  {
    throw Exception("out of range",
    "accessing vector through () - index too small");
  }
  else if (i > mSize)
  {
    throw Exception("length mismatch",
    "accessing vector through () - index too high");
  }
  return mData[i-1];
}

// Member method
// calculate p-norm of a vector v
// default value for p is 2
double Vector::norm(int p) const
{
  double temp, norm_val;

  norm_val = 0.0;
  for (int i=0; i<mSize; i++)
  {
    temp = fabs(mData[i]);
    norm_val += pow(temp, p);
  }
  return pow(norm_val, 1.0/((double) (p)));
}

// Definition of std::cout "<<" operator for a vector
std::ostream& operator<<(std::ostream& output, const Vector& v) 
{
  output << "(";
  for (int i=0; i<v.mSize; i++)
  { 
    output <<  v.mData[i];
    if (i != v.mSize-1)
      output  << ",";
    else
      output  << ")";
  }
  return output;  // for multiple << operators.
}


// return length of a vector
int length(const Vector& v)
{
  return v.mSize;
}

// calculate p-norm of a vector v
// default value for p is 2
double norm(const Vector& v, int p)
{
  double temp, norm_val;

  norm_val = 0.0;
  for (int i=1; i<=length(v); i++)
    {
      temp = fabs(v(i));
      norm_val += pow(temp, p);
    }

  return pow(norm_val, 1.0/((double) (p)));
}

// Fuction that returns as a vector the size of a vector
std::pair<int,int> size(const Vector& v)
{
  return {1,v.mSize};
}
// Vector that cuts a vector of size new_m out of a vector
// if new_m is bigger than the old vector than those entries
// are set to 0.0.
Vector cut(const Vector& v, int new_m)
{
  int min_rows;
  if (v.mSize < new_m)
  {
    min_rows = v.mSize;
  }
  else
  {
    min_rows = new_m;
  }
  Vector new_vec(new_m);
  for (int i = 0; i < min_rows; ++i)
  {
    new_vec.mData[i] = v.mData[i];     
  }
  return new_vec;
}
// Function that writes a vector to the file "voutput.csv"
void writetoCSV(const Vector& v)
{
  std::ofstream out("voutput.csv");
  assert(out.is_open());
  out << v;
  out.close();
  return;
}
