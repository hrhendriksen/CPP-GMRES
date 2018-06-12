#include <iostream>
#include "transport.hpp"
#include <cmath>

Matrix solve(double v, int gridpoints, double endtime, int timepoints, int method)
{
	double delta_x = 1.0/(gridpoints-1);
	double delta_t = endtime/(timepoints-1);

	//-----------------------------------------------------------------
	double sup[gridpoints-1];
	double dia[gridpoints];
	double sub[gridpoints-1];

	// Create S matrix for central difference scheme.
	if(method == 1)
	{
		std::cout << "Solving using central differences \n";
		double a_plus = delta_t/delta_x*(0.5*v-1/delta_x);
		double a_minus = delta_t/delta_x*(-0.5*v-1/delta_x);
		double b = 1+2*delta_t/pow(delta_x,2);


		sup[0]=0;
		dia[0]=1;
		sub[0]=a_plus;

		for (int i = 1; i < gridpoints-2; ++i)
		{
			sup[i] = a_minus;
			dia[i] = b;
			sub[i] = a_plus;
		}

		sup[gridpoints-2] = a_minus;
		dia[gridpoints-2] = b;
		sub[gridpoints-2] = 0;

		dia[gridpoints-1]= 1;
	
	}

	//--------------------------------------------------------------------
	// Create S matrix for upwind difference scheme.
	else if(method == 2)
	{
		std::cout << "Solving using upwind differences \n";
		double d = -delta_t/pow(delta_x,2);
		double e = 1 - v*delta_t/delta_x+2*delta_t/pow(delta_x,2);
		double f = delta_t/delta_x*(v-1/delta_x);

		sup[0]=0;
		dia[0]=1;
		sub[0]=d;

		for (int i = 1; i < gridpoints-2; ++i)
		{
			sup[i] = f;
			dia[i] = e;
			sub[i] = d;
		}

		sup[gridpoints-2] = f;
		dia[gridpoints-2] = e;
		sub[gridpoints-2] = 0;

		dia[gridpoints-1]= 1;
	}
	sparse_trid S(gridpoints,sup,dia,sub);
	print(S);
	//--------------------------------------------------------------------

	// Initial Condition
	Vector u0(gridpoints);
	u0(1)=1;
	
	// Create overall Matrix
	Matrix Z(gridpoints,timepoints);
	// Copy u0 into Z
	Z(1,1)=1;

	double time = 0;
	Vector u_old = u0;

	for (int j = 1; j < timepoints; ++j)
	{
		std::cout << "j is" << j << "\n";
		// Use GMRES to find u1
		// std::cout << "old u_old :"<< u_old << "\n";
		// std::cout << "Steady state??????????" << "\n";
		// std::cout << S_central*u_old-u_old<<"\n";
		Vector u_new = gmres(S, u_old, u_old, gridpoints, 1e-6);
		// std::cout << "u_new :"<< u_new << "\n";

		//Copy u_new into Z
		for (int i = 1; i <= gridpoints; ++i)
		{
			Z(i,j+1) = u_new(i);
		}
		u_old = u_new;
		time += delta_t;
		// std::cout << time << "\n";
	}
	return Z;
}