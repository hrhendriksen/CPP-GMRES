#include <iostream>
#include "gmres.hpp"

Vector gmres(const Matrix& A, Vector& b, Vector& x0, int max_iter, double tol)
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
	Matrix H_hat(max_iter, max_iter+1);
	Matrix V(n,max_iter);
	for (int i = 1; i < n+1; ++i)
	{
		V(i,1) = v(i);
	}
	while(residual > tol && iter < max_iter)
	{		
	 // Do step k of Arnoldi
		Vector w = A * v;
		for (int j = 1; j < iter+1; ++j)
		{
			H_hat(j,iter) = v * w;
			//
			for (int k = 1; k < n+1 ; ++k)
			{
				w(k) -= H_hat(j,iter)*V(k,j);
			}
		}
	H_hat(iter+1,iter) = norm(w);

	for (int i = 1; i < n+1; ++i)
	{
		V(i,iter+1) = w(i)/H_hat(iter+1,iter);
	}

	std::cout << "i is "<<iter<<"\n";
	std::cout << "H_hat is :\n";
	print(H_hat);
	std::cout << "V is :\n";
	print(V);
	iter +=1;
	}
	return v;
}
