#ifndef GIVENSDEF
#define GIVENSDEF

typedef struct {
    double cos, sin;
} angle;

angle Givens(double rho, double sigma);

#endif