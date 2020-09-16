#include <iostream>
#include <tuple>
#include <list>
#include <vector>
#include <functional>
#include "table.h"
#include "makePolynom.h"
#include "function.h"
using namespace std;

int main()
{
    cout << "Polynomial interpolation problem:\n";
    cout << "Variant 13:\n";

    int m = 0;
    double a = 0;
    double b = 0;
    double x = 0;
    int n = 0;

    cout << "Enter number of values in table (m + 1):\n";
    cin >> m;
    --m;
    cout << "m:" << m << "\n";

    cout << "Enter a:\n";
    cin >> a;

    cout << "Enter b:\n";
    cin >> b;

    printf("[a, b]: [%.14lf, %.14lf]", a, b);

    char key;
    list<tuple<double, double>> table = makeTable(a, b, m);
    do
    {
        cout << "Enter interplation polynom degree n, n <= " << m << "\n";
        cin >> n;
        while (n > m)
        {
            cout << "Entered an invalid value of n, n should be less than " << m << ", try again :)" << "\n";
            cin >> n;
        }

        cout << "Enter x:\n";
        cin >> x;
        printf("x: %.14lf\n", x);

        printf("[Xk, f(Xk)]\n");
        list<tuple<double, double>> sortedTable = sortTable(table, x);

        double result = polynom(sortedTable, n, x);
        printf("%.14lf \n", result);

        double prec = precision(result, x);
        printf("precision: %.14lf \n", prec);

        cout << "Enter \"q\" to quit, enter any key to continue with new n and x:\n";
        cin >> key;
    } while (key != 'q');

    return 0;
}