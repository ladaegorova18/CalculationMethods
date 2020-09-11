#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cmath>
#include <list>
#include <tuple>
#include "RootsSeparation.h"
#include "Bisection.h"
#include "NewtonMethod.h"
#include "modifiedNewtonMethod.h"
#include "secantMethod.h"

const double epsilon = 0.00000000001;

int main()
{
    double A = 0;
    double B = 0;
    int N = 0;

    std::cout << "Enter A:\n";
    std::cin >> A;
    std::cout << "Enter B:\n";
    std::cin >> B;
    std::cout << "Enter N:\n";
    std::cin >> N;

    std::cout << "Roots separation:\n";
    list<tuple<double, double>> intervals = separation(A, B, N);

    printf("Epsilon: %.10lf\n", epsilon);

    for (auto inter : intervals)
    {

        double newA = 0;
        double newB = 0;
        tie(newA, newB) = inter;
        printf("Roots in interval [%f, %f]\n", newA, newB);

        bisection(newA, newB, epsilon);
        std::cout << "\n";

        newtonMethod(newA, newB, epsilon);
        std::cout << "\n";

        modifiedNewtonMethod(newA, newB, epsilon);
        std::cout << "\n";

        secantMethod(newA, newB, epsilon);
        std::cout << "\n";
    }
    return 0;
}