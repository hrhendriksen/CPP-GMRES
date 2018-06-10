#ifndef GMRESDEF
#define GMRESDEF

#include <cmath>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

typedef struct {
    double cos, sin;
} angle;

Vector gmres(const Matrix& A, Vector& b, Vector& x0, int max_it, double tol);

angle Givens(double rho, double sigma);


#endif
