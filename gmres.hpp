#ifndef GMRESDEF
#define GMRESDEF

#include <cmath>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

typedef struct {
    double cos, sin;
} angle;

angle Givens(double rho, double sigma);

 
template <class matrix_type>
Vector gmres(const matrix_type& A, const Vector& b, const Vector& x0, int max_it, double tol);

#include "gmres.tcc"

#endif
