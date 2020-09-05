#include "NewtonMethod.h"
#include "function.h"
#include <cmath>

void secantMethod(double A, double B, double epsilon)
{
	printf("Secant method:\n");
	int N = 0;
	double Xk = A;
	double Xk_1 = B;
	double Xk_2 = 0;
	double X0 = Xk - (function(Xk) * (Xk_1 - Xk) / (function(Xk_1) - function(Xk)));

	do
	{
		Xk_2 = Xk - (function(Xk) * (Xk_1 - Xk) / (function(Xk_1) - function(Xk)));
		Xk = Xk_1;
		Xk_1 = Xk_2;
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