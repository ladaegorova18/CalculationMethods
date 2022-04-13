#include "function.h"

double function(double x)
{
    return x - cos(M_PI * x) * cos(M_PI * x);
}

double fstDerivative(double x)
{
    return 1 + 2 * M_PI * sin(M_PI * x) * cos(M_PI * x);
}

double sndDerivative(double x)
{
    return 2 * M_PI * M_PI * cos(2 * M_PI * x);
}
