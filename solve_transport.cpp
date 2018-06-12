#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "gmres.hpp"
#include "transport.hpp"
#include <time.h>       /* time */
#include <typeinfo>

int main(int argc, char const *argv[])
{
	//print(solve(-10,100,1,100,1)); Good agreement with exact solution.
	//print(solve(-10,200,1,200,1));
	// print(solve(-20,100,1,100,1));
	// Matrix A(3,3);
	// std::cout << A;
	writetoCSV(solve(-100,100,1,100,1));
	return 0;
}