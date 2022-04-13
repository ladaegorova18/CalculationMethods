#include "RootsSeparation.h"
#include "function.h"
#include <inttypes.h>

list<tuple<double, double>> separation(double A, double B, int N)
{
    list<tuple<double, double>> intervals;
    double H = (B - A) / N;
    int counter = 0;
    double x1 = A;
    double x2 = x1 + H;
    double y1 = function(x1);
    double y2 = 0;

    while (x2 <= B)
    {
        y2 = function(x2);
        if (y1 * y2 <= 0)
        {
            ++counter;
            tuple<double, double> interval(x1, x2);
            intervals.push_back(interval);
            printf("[%f, %f]\n", x1, x2);
        }
        x1 = x2;
        x2 = x1 + H;
        y1 = y2;
    }
    cout << counter << "\n";
    return intervals;
}
