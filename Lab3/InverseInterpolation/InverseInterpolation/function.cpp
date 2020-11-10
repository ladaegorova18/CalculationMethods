#include "function.h"

double function(double x)
{
	return log(1 + x) - exp(x);
}

double invPrecision(double result, double x)
{
	return abs(function(result) - x);
}