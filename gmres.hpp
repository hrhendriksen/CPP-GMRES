#ifndef GMRESDEF
#define GMRESDEF

#include <cmath>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

typedef struct {
    double cos, sin;
} angle;
 
template <class matrix_type>
Vector gmres(const matrix_type& A, Vector& b, Vector& x0, int max_it, double tol);
// Vector gmres(const sparse_trid& A, Vector& b, Vector& x0, int max_iter, double tol);

angle Givens(double rho, double sigma);


#endif
