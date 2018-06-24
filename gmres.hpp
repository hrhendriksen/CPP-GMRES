#ifndef GMRESDEF
#define GMRESDEF

#include <cmath>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Givens.hpp"

template <class matrix_type>
Vector gmres(const matrix_type& A, const Vector& b, const Vector& x0, int max_it, double tol);

#include "gmres.tcc"
#endif
