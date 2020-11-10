#include "PolynomialThirdDegree.h"
#include <cmath>

double PolynomialThirdDegree::func(double x)
{
    return pow(x, 3) * (-17) + 0.3 * pow(x, 2) - x + 5;
}

double PolynomialThirdDegree::integral(double x)
{
    return pow(x, 4) * (-4.25) + pow(x, 3) * 0.1  - pow(x, 2) / 2 + 5 * x;
}

double PolynomialThirdDegree::derivative(double x)
{
    return -54 * pow(x, 2) + 0.6 * x - 1;
}

double PolynomialThirdDegree::sndDerivative(double x)
{
    return -108 * x + 0.6;
}

double PolynomialThirdDegree::fourthDerivative(double x)
{
    return 0.0;
}