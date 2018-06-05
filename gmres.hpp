#ifndef GMRESDEF
#define GMRESDEF

#include <cmath>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

Vector gmres(const Matrix& A, Vector& b, Vector& x0, int max_it, double tol);

#endif
