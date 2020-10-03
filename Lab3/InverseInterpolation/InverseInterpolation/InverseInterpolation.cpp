#include <iostream>
#include <list>
#include "table.h"
#include "newton.h"
#include "function.h"
#include "InverseInterpolation.h"
using namespace std;

int main()
{
    cout << "Inverse interpolation problem:\n";
    cout << "Variant 13:\n";

    double a = 0;
    double b = 1;
    double x = 0;
    double epsilon = 0.00000001;
    double F = 0;
    int n = 0;
    int m = 10;

    char key;
    list<tuple<double, double>> table = makeTable(a, b, m);
    do
    {
        cout << "Enter F:\n";
        cin >> F;
        printf("F: %.14lf\n", F);

        cout << "Enter interplation polynom degree n, n <= " << m << "\n";
        cin >> n;
        while (n > m)
        {
            cout << "Entered an invalid value of n, n can not be greater than " << m << ", try again :)" << "\n";
            cin >> n;
        }

        printf("Epsilon is %.14lf, do you want to change epsilon? (y/n)\n", epsilon);
        cin >> key;
        if (key == 'y')
        {
            cout << "Enter epsilon:\n";
            cin >> epsilon;
            printf("epsilon: %.14lf\n", epsilon);
        }

        firstMethod(table, F, n);
        secondMethod(table, F, n, a, b, epsilon);

        cout << "Enter \"q\" to quit, enter any key to continue with new F, n and epsilon:\n";
        cin >> key;
    } while (key != 'q');

    return 0;
}

void firstMethod(std::list<std::tuple<double, double>>& table, double F, int n)
{
    cout << "First method:\n";

    list<tuple<double, double>> swappedTable = swapTable(table);

    cout << "Sorted table:\n";
    printf("Xk \t\t\t f(Xk)\n");
    list<tuple<double, double>> sortedTable = sortTable(swappedTable, F);
    auto newton = new Newton(sortedTable, n);

    double result = newton -> count(F);
    printf("Function value in x: %.14lf \n", result);

    double prec = invPrecision(result, F);
    printf("precision: %.14lf \n", prec);
}

void secondMethod(std::list<std::tuple<double, double>>& table, double F, int n, double a, double b, double epsilon)
{
    cout << "Second method:\n";

    list<tuple<double, double>> sortedTable = modifiedSort(table, F);
    auto newton = new Newton(sortedTable, n);
    double result = newton -> bisection(a, b, epsilon, F);
    printf("Function value in x: %.14lf \n", result);

    double prec = invPrecision(result, F);
    printf("precision: %.14lf \n", prec);
}
