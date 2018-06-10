#include <iostream>
#include "gmres.hpp"
#include <cmath>

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

	// Vector v_iter(n);
	Matrix H(max_iter, max_iter+1);
	Matrix V(n,max_iter);
	for (int i = 1; i < n+1; ++i)
	{
		V(i,1) = v(i);
	}
	//Create a vector for the cosine and sine
	Vector cosine(max_iter);
	Vector sine(max_iter);

	while(residual > tol && iter < max_iter)
	{		
	 // Do step k of Arnoldi

		Vector v_iter(n);
		for (int i = 1; i < n+1; ++i)
		{
			v_iter(i) = V(i,iter);
		}
		Vector w = A * v_iter;

		for (int j = 1; j < iter+1; ++j)
		{
			H(j,iter) = v * w;
			//
			for (int k = 1; k < n+1 ; ++k)
			{
				w(k) -= H(j,iter)*V(k,j);
			}
		}
	H(iter+1,iter) = norm(w);

	for (int i = 1; i < n+1; ++i)
	{
		V(i,iter+1) = w(i)/H(iter+1,iter);
	}

	std::cout << "i is "<<iter<<"\n";
	std::cout << "H is :\n";
	print(H);
	// std::cout << "V is :\n";
	// print(V);
	std::cout << "Now solve LLS:\n";

	// Apply Gives rotation to H
	double temp;
	for (int i = 1; i < iter; ++i)
	{
		temp = cos(i)*H(i,iter)+sin(i)*H(i+1,iter);
		H(i+1, iter) = -sin(i) * H(i,iter)+cos(i)*H(i+1,iter);
		H(i,iter) = temp;
	}

	double rho;
	rho = H(iter,iter);

	double sigma;
	sigma = H(iter+1,iter);

	angle ang  = Givens(rho, sigma);
	// std::cout<<"Rho, Sigma"<< rho<< " , "<< sigma<<"\n";
	// std::cout<<"Hypothetical angles"<<ang.cos << " and " <<ang.sin<<"\n";
	cosine(iter) = ang.cos;
	sine(iter) = ang.sin;


	// Recent Givens rotation
	H(iter,iter) = cos(iter)*H(iter,iter)+sin(iter)*H(iter+1,iter);
	// Now we must have a triangular H
	std::cout << "Now H triangular \n";
	print(H);
	iter +=1;
	}
	return x0;
}

angle Givens(double rho, double sigma) {
	double cos;
	double sin;
	if (rho == 0)
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


