#include <iostream>
#include "gmres.hpp"

Vector gmres(const Matrix& A, Vector& b, Vector& x0, int max_it, double tol)
{
	assert(A.GetNumberofRows() == A.GetNumberofColumns());

	int iter = 0;

	// Calculate initial residual r0
	Vector r0 = b - A*x0;
	double norm_r0 = norm(r0);

	if (norm(r0)<tol)
	{
		return x0;
	}

	Vector v = r0/norm_r0;

	//Now build the matrix V using the previously found v's.
	Matrix V(v);
	print(V);
	return x0;
}