#include <iostream>
#include <cmath>

template <class matrix_type>
Vector gmres(const matrix_type& A, Vector& b, Vector& x0, int max_iter, double tol)
{
	// Assert that A is a square matrix
    int n = A.GetNumberofRows();
	assert(n == A.GetNumberofColumns());

	// Initialise an iteration counter
	int iter = 1;

	// Preallocate and (possibly) initialise all used vectors and matrices:
	// Calculate initial residual r0
	Vector r0 = b - A*x0;
	double norm_r0 = norm(r0);
	double norm_b = norm(b);
	double error = norm_r0/norm_b;

	// Preallocate an overall vector for the cosine and sine of the Givens rotations
	Vector cosine(max_iter);		
	Vector sine(max_iter);

	// Preallocate and initialise beta vector, beta = (r_0,0,0,0,0,0,0,0)
	Vector beta(n+1);
	beta(1) = norm_r0;

	// Preallocate and initialise an overall y and residuals vector
	Vector residuals(max_iter);
	residuals(1) = error;
	Vector y(max_iter);


	if (error<tol)
	{
		return x0;
	}

	// Preallocate an overall H and V matrix, set first column of V to v1
	Matrix H(max_iter, max_iter+1);
	Matrix V(n,max_iter);

	Vector v_1 = r0/norm_r0;

	for (int i = 1; i < n+1; ++i)
	{
		V(i,1) = v_1(i);
	}

	// Start the GMRES iteration until the desired tolerance is reached 
	// or the maximum number of iterations have been executed.

	while(error > tol && iter < max_iter)
	{		
		//Take the iter^{th} column of V
		Vector v_iter(n);

		for (int i = 1; i < n+1; ++i)
		{
			v_iter(i) = V(i,iter);
		}

		// Do step k of Arnoldi
		Vector w = A * v_iter;

	// 	// Fill in the iter^{th} column of H
		for (int j = 1; j < iter+1; ++j)
		{
			for (int i = 1; i < n+1; ++i)
			{
				H(j,iter) += w(i)*V(i,j);
			}
			
			for (int k = 1; k < n+1 ; ++k)
			{
				w(k) -= H(j,iter)*V(k,j);
			}
		}
		H(iter+1,iter) = norm(w);

		// Fill in the iter^{th} column of V
		for (int i = 1; i < n+1; ++i)
		{
			V(i,iter+1) = w(i)/H(iter+1,iter);
		}

		// Apply Gives rotation to H
		double temp;

		// Firstly, work with the (iter-1) previous Givens rotations on the
		// upper (iter) entries of the iter^{th} column of H
		for (int i = 1; i < iter; ++i)
		{
			temp = cosine(i)*H(i,iter)+sine(i)*H(i+1,iter);
			H(i+1, iter) = -sine(i) * H(i,iter)+cosine(i)*H(i+1,iter);
			H(i,iter) = temp;
		}

		// Now find the Givens rotation for the (iter,iter) and (iter+1, iter)
		// entries of H.
		double rho;
		rho = H(iter,iter);

		double sigma;
		sigma = H(iter+1,iter);

		// Use the function Givens to calculate the cos and sin and returns a struc
		angle ang  = Givens(rho, sigma);

		// Add these angles to your cosine and sine vector
		cosine(iter) = ang.cos;
		sine(iter) = ang.sin;

		// Apply the recent Givens rotation
		H(iter,iter) = cosine(iter)*H(iter,iter)+sine(iter)*H(iter+1,iter);
		H(iter+1,iter) = 0;

		// Also work with our Givens rotation on the beta = ||r_0||*e_1 vector
		beta(iter+1) = -sine(iter)*beta(iter);
		beta(iter) = cosine(iter)*beta(iter);

		//The error will be the absolute value of the last entry of our beta vector
		error = std::fabs(beta(iter+1))/norm_b;
	 	residuals(iter+1) = error;
	 
	 	// Now copy the upper (iter) entries of beta; g_n
	 	// to apply backward substitution.
		Vector g_n = cut(beta, iter);	

	// 	// Backward substitution of the triangular system: R*y = g_n
		for (int i=iter; i>=1; i--)
		{
			y(i) = g_n(i)/H(i,i);

			for (int k=i-1; k>=1; k--)
			{	
				g_n(k) -= H(k,i)*y(i);
			}
		}

		iter +=1;
	}

	// Calculate the final solution
	Vector delta_x(n);

	 for (int i = 1; i <= n; ++i)
	 {
	 	for (int j = 1; j <= iter-1; ++j)
	 	{
	 		delta_x(i) += V(i,j)*y(j);
	 	}
	 }
	// std::cout<<cut(residuals,iter)<<"\n";

	// Choose to either let the function return the solution or the vector of residuals
	return x0 + delta_x;
	 // return residuals;
}

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