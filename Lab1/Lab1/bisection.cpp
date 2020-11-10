#include "Bisection.h"
#include "function.h"
#include <cmath>

void bisection(double a, double b, double epsilon)
{
	printf("Bisection method:\n");
	int N = 0;
	double X0 = (a + b) / 2;

	do
	{
		double c = (a + b) / 2;
		if (function(a) * function(c) <= 0)
		{
			b = c;
		}
		else
		{
			a = c;
		}
		++N;
	} 
	while (b - a > 2 * epsilon);
	
	double X = (a + b) / 2;
	double delta = (b - a) / 2;
	
	printf("X0: %.14lf\n", X0);
	printf("Delta: %.14lf\n", delta);
	printf("Steps: %d\n", N);
	printf("X: %.14lf\n", X);
	printf("|Xn - Xn-1:|: %.14lf\n", abs(b - a));
	printf("|f(X) - 0|: %.14lf\n", abs(function(X) - 0.0));
}
