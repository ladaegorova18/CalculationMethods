#define _USE_MATH_DEFINES
#include <cmath>
#include "Function.h"
#include <functional>
using namespace std;

double Function::p(double x)
{
    return pow(M_E, -x);
}

double Function::f(double x)
{
    return sin(x);
}

double Function::func(double x)
{
    return sin(x) * pow(M_E, -x);
}

double Function::p_1(double x)
{
    return pow(M_E, -x) * x;
}

double Function::p_2(double x)
{
    return pow(M_E, -x) * pow(x, 2);
}

double Function::p_3(double x)
{
    return pow(M_E, -x) * pow(x, 3);
}

double Function::integral(double x)
{
    return -0.5 * pow(M_E, -x) * (sin(x) + cos(x));
}

double Function::precision(double result1, double result2)
{
    return abs(result1 - result2);
}


//double function(double) Function::multiplyF(double x, int i)
//{
//    double function(double x)
//    {
//        sin(x)* pow(x, i)
//    };
//    return function;
//}
