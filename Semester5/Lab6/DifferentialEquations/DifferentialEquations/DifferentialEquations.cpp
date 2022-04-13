#include <iostream>
#include <clocale>
#include <iomanip>
#include <Windows.h>
#include "Solver.h"
#include "DifferentialEquations.h"
using namespace std;

void printFirstLine()
{
    printOneString("x");
    printOneString("y");
    printOneString("Погрешность");
    cout << "|" << endl;
}

void printOneString(string letter)
{
    cout << "|" << setw(20) << letter;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Численное решение Задачи Коши для обыкновенного дифференциального уравнения первого порядка\n";
    cout << "Вариант 13:\n";
    cout << "y' = -y + y^2 + 1\n";

    double h;
    int N;
    double y0 = 0.25; 
    double x0 = 0;
    vector<double> X;
    vector<double> Y;
    char key = ' ';

    do
    {
        cout << "Введите интервал деления h: ";
        cin >> h;
        cout << "\n";

        cout << "Введите число точек правее от х0 (N): ";
        cin >> N;
        cout << "\n";

        Solver* solver = new Solver(0, y0, h, N);

        cout << "\x1b[32mМетод Тейлора:\x1b[0m" << endl;
        for (int k = -2; k <= N; ++k)
        {
            double x = x0 + k * h;
            double preciseSolution = solver->presiseSolution(x);

            X.push_back(x);
            Y.push_back(preciseSolution);

            double y = solver->Taylor(x);
            printValue(X[k + 2], y, Y[k + 2], solver);
            if (k == N) printYn(y, Y[k + 2]);
        }
        cout << endl;

        cout << "\x1b[32mМетод Адамса 4-го порядка:\x1b[0m" << endl;

        vector<double> yAdams = solver->Adams(N);
        for (int k = 3; k <= N; ++k)
        {
            printValue(X[k + 2], yAdams[k + 2], Y[k + 2], solver);
            if (k == N) printYn(yAdams[k + 2], Y[k + 2]);
        }
        cout << endl;

        cout << "\x1b[32mМетод Рунге-Кутта:\x1b[0m" << endl;
        for (int k = 1; k <= N; ++k)
        {
            double yRungeKutta = solver->RungeKutta(k);
            printValue(X[k + 2], yRungeKutta, Y[k + 2], solver);
            if (k == N) printYn(yRungeKutta, Y[k + 2]);
        }
        cout << endl;

        cout << "\x1b[32mМетод Эйлера:\x1b[0m" << endl;
        for (int k = 1; k <= N; ++k)
        {
            double yEuler = solver->Euler(X[k + 2], k);
            printValue(X[k + 2], yEuler, Y[k + 2], solver);
            if (k == N) printYn(yEuler, Y[k + 2]);
        }
        cout << endl;

        cout << "\x1b[32mМетод Эйлера I:\x1b[0m" << endl;
        for (int k = 1; k <= N; ++k)
        {
            double yEulerI = solver->EulerI(X[k + 2], k);
            printValue(X[k + 2], yEulerI, Y[k + 2], solver);
            if (k == N) printYn(yEulerI, Y[k + 2]);
        }
        cout << endl;

        cout << "\x1b[32mМетод Эйлера II:\x1b[0m" << endl;
        for (int k = 1; k <= N; ++k)
        {
            double yEulerII = solver->EulerII(X[k + 2], k);
            printValue(X[k + 2], yEulerII, Y[k + 2], solver);
            if (k == N) printYn(yEulerII, Y[k + 2]);
        }
        cout << endl;

        cout << "Введите q, чтобы выйти:\n";
        cin >> key;
    } while (key != 'q');
}

void printYn(double y, double Yn)
{
    cout << setprecision(14) << "   |y(xN) - yN|: " << fabs(y - Yn) << endl;
}

void printValue(double x, double y, double exactValue, Solver* solver)
{
    double precision = solver->precision(y, exactValue);
    cout << setprecision(14) << "|" << setw(20) << x << "|" << setw(20) << y << "|" << setw(20) << precision << "|" << "\n";
}

