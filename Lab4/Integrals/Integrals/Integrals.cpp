#include <iostream>
#include "Integral.h"
#include <locale>
#include "PolynomialFirstDgr.h"
#include "PolynomialThirdDegree.h"
#include "PolynomialZeroDgr.h"
#include "Exponenta.h"
#include "SomeFunction.h"
using namespace std;

int main()
{
    setlocale(LC_ALL, "Russian");
    double A = 0;
    double B = 0;
    int m = 0;
    double h = 0;

    cout << "Приближённое вычисление интеграла по составным квадратурным формулам \n";

    char key = 'w';
    do
    {
        cout << "f(x) = -5.5\n";
        cout << "f(x) = 2.3 * x - 3\n";
        cout << "f(x) = x ^ 3 * (-17) + 0.3 * x ^ 2 - x + 5\n";
        cout << "f(x) = e^x\n";
        cout << "f(x) = 5 * cos(x) + 3\n";

        cout << "Введите A:\n";
        cin >> A;
        cout << "Введите B:\n";
        cin >> B;
        cout << "Введите число промежутков деления m:\n";
        cin >> m;
        h = (B - A) / m;

        cout << "A: " << A << "\n";
        cout << "B: " << B << "\n";
        cout << "m: " << m << "\n";
        cout << "h: " << h << "\n";

        cout << "f(x) = -5.5\n";
        PolynomialZeroDgr* funZeroDgr = new PolynomialZeroDgr();
        Integral* integralZeroDgr = new Integral(*funZeroDgr, A, B, m);
        integralZeroDgr->countIntegral();

        cout << "\n";
        cout << "f(x) = 2.3 * x - 3\n";
        PolynomialFirstDgr* funFirstDgr = new PolynomialFirstDgr();
        Integral* integralFirstDgr = new Integral(*funFirstDgr, A, B, m);
        integralFirstDgr->countIntegral();

        cout << "\n";
        cout << "f(x) = x ^ 3 * (-17) + 0.3 * x ^ 2 - x + 5\n";
        PolynomialThirdDegree* funThirdDgr = new PolynomialThirdDegree();
        Integral* integralThirsDgr = new Integral(*funThirdDgr, A, B, m);
        integralThirsDgr->countIntegral();

        cout << "\n";
        cout << "f(x) = e^x\n";
        Exponenta* exp = new Exponenta();
        Integral* integralExp = new Integral(*exp, A, B, m);
        integralExp->countIntegralWithTeorPrecision();

        cout << "\n";
        cout << "f(x) = 5 * cos(x) + 3\n";
        SomeFunction* someFunc = new SomeFunction();
        Integral* integralSomeFunc = new Integral(*exp, A, B, m);
        integralSomeFunc->countIntegral();

        cout << "Нажмите q, чтобы выйти:\n";
        cin >> key;

    } while (key != 'q');
}
