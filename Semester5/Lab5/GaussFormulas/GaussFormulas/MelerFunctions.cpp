#define _USE_MATH_DEFINES
#include <cmath>
#include "MelerFunctions.h"

double MehlerFunctions::p(double x)
{
    return 1.0 / (sqrt(1.0 - pow(x, 2)));
}

double MehlerFunctions::f(double x)
{
    return sin(x);
}
