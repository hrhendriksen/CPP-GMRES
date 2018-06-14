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
// Helper function to calculate the entries of the Givens matrix
angle Givens(double rho, double sigma) {
	double cos;
	double sin;
	if (rho == 0.0)
	{
		cos = 0;
		sin = 1;
	}
	else
	{
		double total = sqrt(pow(rho,2)+pow(sigma,2));
		cos = std::fabs(rho)/total;
		sin = cos*sigma/rho;
	}

    angle stru = {cos, sin};
    return stru;
}
 
template <class matrix_type>
Vector gmres(const matrix_type& A, Vector& b, Vector& x0, int max_it, double tol);

#include "gmres.tcc"

#endif
