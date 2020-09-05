#include "NewtonMethod.h"
#include "function.h"
#include <cmath>

void newtonMethod(double A, double B, double epsilon)
{
	printf("Newton method:\n");
	int N = 0;
	double Xk = B;
	double Xk_1 = 0;
	double X0 = 0;

	if (function(A) * sndDerivative(A) > 0)
	{
		X0 = A;
	}
	else if (function(B) * sndDerivative(B) > 0)
	{
		X0 = B;
	}
	else X0 = (A + B) / 2;

	Xk = X0;
	do
	{
		Xk = Xk_1;
		Xk_1 = Xk - function(Xk) / fstDerivative(Xk);
		++N;
	} while (abs(Xk_1 - Xk) > epsilon);

	double X = (Xk + Xk_1) / 2;
	double delta = (abs(Xk_1 - Xk)) / 2;

	printf("Delta: %.10lf\n", delta);
	printf("X0: %.10lf\n", X0);
	printf("Steps: %d\n", N);
	printf("X: %.10lf\n", X);
	printf("|Xn - Xn-1:|: %.10lf\n", abs(Xk_1 - Xk));
	printf("|f(X) - 0|: %.10lf\n", abs(function(X) - 0.0));
}