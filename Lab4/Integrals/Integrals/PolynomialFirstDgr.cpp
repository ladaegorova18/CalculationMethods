#include "PolynomialFirstDgr.h"
#include <cmath>

double PolynomialFirstDgr::func(double x)
{
    return 2.3 * x - 3;
}

double PolynomialFirstDgr::integral(double x)
{
    return 1.15 * pow(x, 2) - 3 * x;
}

double PolynomialFirstDgr::derivative(double x)
{
    return 2.3 * x;
}

double PolynomialFirstDgr::sndDerivative(double x)
{
    return 0.0;
}

double PolynomialFirstDgr::fourthDerivative(double x)
{
    return 0.0;
}
