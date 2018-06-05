#include <iostream>
#include "gmres.hpp"

Vector gmres(const Matrix& A, Vector& b, Vector& x0, int max_it, double tol)
{
    int n = A.GetNumberofRows();
	assert(n == A.GetNumberofColumns());

	int iter = 1;
	double residual = 5.0;
	// Calculate initial residual r0
	Vector r0 = b - A*x0;
	double norm_r0 = norm(r0);

	if (norm(r0)<tol)
	{
		return x0;
	}

	Vector v = r0/norm_r0;

	Vector v_iter(n);

	while(residual > tol && iter < max_it)
	{
		Matrix H(iter, iter+1);
	 // Do step k of Arnoldi
		Vector w = A * v;
		for (int j = 1; j < iter; ++j)
		{
			continue;
		}




		iter +=1;
	}
	return v;
}
