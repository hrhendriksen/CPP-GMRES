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
	// Function that solves our transport problem with Dirichlet BC
	// and directly writes it into the output.CSV file.
	// The last argument lets the user choose between 3 options;
	
	// 1: Central Differences Sparse Matrix
	// 2: Upwind Differences Sparse Matrix
	// 3: Dense Upwind Differences Matrix

	writetoCSV(solve(0,100,1,100,2));

	return 0;
}