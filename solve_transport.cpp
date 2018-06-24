#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "transport.hpp"

int main(int argc, char const *argv[])
{
	// Function that solves our transport problem with Dirichlet BC
	// and directly writes it into the output.CSV file.
	// The last argument lets the user choose between 2 options;
	
	// 1: Central Differences Sparse Matrix
	// 2: Upwind Differences Sparse Matrix
	writetoCSV(solve(100,100,1,100,2));
	return 0;
}