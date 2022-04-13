#include "NewtonMethod.h"
#include "function.h"
#include <cmath>

void modifiedNewtonMethod(double A, double B, double epsilon)
{
	printf("Modified Newton method:\n");

	int N = 0;
	double X0 = (A + B) / 2;
	double X1 = A;
	double X2 = X0;

	double fstDrvX0 = fstDerivative(X0);
	while (abs(X2 - X1) > epsilon)
	{
		X1 = X2;
		X2 = X1 - function(X1) / fstDrvX0;
		++N;
	}

	double X = (X1 + X2) / 2;

	double delta = (abs(X2 - X1)) / 2;

	printf("Delta: %.14lf\n", delta);
	printf("X0: %.14lf\n", X0);
	printf("Steps: %d\n", N);
	printf("X: %.14lf\n", X);
	printf("|Xn - Xn-1:|: %.14lf\n", abs(X2 - X1));
	printf("|f(X) - 0|: %.14lf\n", abs(function(X) - 0.0));
}
