#include "Givens.hpp"
#include <cmath>

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