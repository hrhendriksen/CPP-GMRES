/*
 * Q9.hpp
 *
 *  Created on: 25 apr. 2018
 *      Author: hiddehendriksen
 */
#include "Matrix.hpp"

Matrix& create_aug(const Vector& v, const Matrix& m);
double** Gaussian_elimination(double** aug, int n);
void solve_triangular(double** aug, int n);
double** allocatematrix(int size_a, int size_b);
class Matrix;