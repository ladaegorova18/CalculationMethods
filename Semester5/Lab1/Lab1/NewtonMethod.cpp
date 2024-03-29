#include "NewtonMethod.h"
#include "function.h"
#include <cmath>

void newtonMethod(double A, double B, double epsilon)
{
	printf("Newton method:\n");
	int N = 0;
	double X1 = 0;
	double X0 = (A + B) / 2;
	double X2 = X0;

	do
	{
		int p = 1;
		X1 = X2;
		if (fstDerivative(X1) == 0)
		{
			p += 2;
			X2 = X0;
		}
		else
		{
			X2 = X1 - p * function(X1) / fstDerivative(X1);
		}
		++N;
	} while (abs(X2 - X1) > epsilon);

	double X = (X1 + X2) / 2;

	double delta = (abs(X2 - X1)) / 2;

	printf("Delta: %.14lf\n", delta);
	printf("X0: %.14lf\n", X0);
	printf("Steps: %d\n", N);
	printf("X: %.14lf\n", X);
	printf("|Xn - Xn-1:|: %.14lf\n", abs(X2 - X1));
	printf("|f(X) - 0|: %.14lf\n", abs(function(X) - 0.0));
}