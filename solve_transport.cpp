#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "Exception.hpp"
#include "Matrix.hpp"
// #include "gmres.hpp"
#include "transport.hpp"

int main(int argc, char const *argv[])
{
	// Function that solves our transport problem with Dirichlet BC
	// and directly writes it into the output.CSV file.
	// The last argument lets the user choose between 3 options;
	
	// 1: Central Differences Sparse Matrix
	// 2: Upwind Differences Sparse Matrix
	// 3: Dense Upwind Differences Matrix
	
	// clock_t t;
	// t = clock();
	writetoCSV(solve(0,100,1,100,2));
	// t = clock() - t;
	// printf ("It took me %lu clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	std::cout << "asdf";
	return 0;
}