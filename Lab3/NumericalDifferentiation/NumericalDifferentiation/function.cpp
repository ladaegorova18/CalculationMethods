#include "function.h"

double function(double x)
{
	return exp(6 * x);
}

double fstDerivative(double x)
{
	return 6 * exp(6 * x);
}

double sndDerivative(double x)
{
	return 36 * exp(6 * x);
}

double precision(double result, double x, double func(double))
{
	return abs(result - func(x));
}