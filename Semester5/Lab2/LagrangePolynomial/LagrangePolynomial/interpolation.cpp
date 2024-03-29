#include <iostream>
#include <tuple>
#include <list>
#include <vector>
#include <functional>
#include "table.h"
#include "lagrangePolynomial.h"
#include "newtonPolynomial.h"
#include "function.h"
using namespace std;

int main()
{
    cout << "Polynomial interpolation problem:\n";
    cout << "Variant 13:\n";

    int m_plus_1 = 0;
    double a = 0;
    double b = 0;
    double x = 0;
    int n = 0;

    cout << "Enter number of values in table (m + 1):\n";
    cin >> m_plus_1;
    int m = m_plus_1 - 1;

    cout << "m:" << m << "\n";

    cout << "Enter a:\n";
    cin >> a;

    cout << "Enter b:\n";
    cin >> b;
    while (a >= b)
    {
        cout << "Entered an invalid value of b, b should be greater than " << a << ", try again :)" << "\n";
        cin >> b;
    }

    printf("[a, b]: [%.14lf, %.14lf]", a, b);
    cout << "\n";

    char key;
    list<tuple<double, double>> table = makeTable(a, b, m);
    do
    {
        cout << "Enter interplation polynom degree n, n <= " << m << "\n";
        cin >> n;
        while (n > m)
        {
            cout << "Entered an invalid value of n, n can not be greater than " << m << ", try again :)" << "\n";
            cin >> n;
        }

        cout << "Enter x:\n";
        cin >> x;
        printf("x: %.14lf\n", x);

        cout << "Sorted table:\n";
        printf("Xk \t\t\t f(Xk)\n");
        list<tuple<double, double>> sortedTable = sortTable(table, x);

        cout << "Lagrange method:\n";
        double lagrangeResult = polynomial(sortedTable, n, x);
        printf("Function value in x: %.14lf \n", lagrangeResult);

        double prec = precision(lagrangeResult, x);
        printf("precision: %.14lf \n", prec);

        cout << "Newton method:\n";
        auto newtonClass = new Newton(sortedTable, n);
        double newtonResult = newtonClass -> count(x);
        printf("Function value in x: %.14lf \n", newtonResult);

        double prec1 = precision(newtonResult, x);
        printf("precision: %.14lf \n", prec1);

        cout << "Enter \"q\" to quit, enter any key to continue with new n and x:\n";
        cin >> key;
    } while (key != 'q');

    return 0;
}