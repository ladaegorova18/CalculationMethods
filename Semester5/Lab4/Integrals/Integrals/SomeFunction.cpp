#include "SomeFunction.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include<iostream>
#include<iomanip>  

double SomeFunction::func(double x)
{
    return 5 * cos(x) + 3;
}

double SomeFunction::integral(double x)
{
    return 5 * sin(x) + 3 * x;
}

double SomeFunction::derivative(double x)
{
    return -5 * sin(x);
}

double SomeFunction::sndDerivative(double x)
{
    return -5 * cos(x);
}

double SomeFunction::fourthDerivative(double x)
{
    return 5 * cos(x);
}