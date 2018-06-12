#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "Exception.hpp"
#include "Matrix.hpp"
#include "gmres.hpp"
#include "transport.hpp"
#include <time.h>       /* time */

int main(int argc, char const *argv[])
{
	Matrix A = solve(10,10,1,10);
	print(A);
	return 0;
}