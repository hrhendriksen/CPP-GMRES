#ifndef GMRESDEF
#define GMRESDEF

#include <cmath>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Exception.hpp"

Vector gmres(const Matrix& A, Vector& b, Vector& x0, int max_it, double tol);

#endif