#ifndef TRANSPORTDEF
#define TRANSPORTDEF

#include <cmath>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "gmres.hpp"
#include <iostream>

Matrix solve(double v=2, int gridpoints=10, double endtime=10, int timepoints=10, int method = 1);

#endif
