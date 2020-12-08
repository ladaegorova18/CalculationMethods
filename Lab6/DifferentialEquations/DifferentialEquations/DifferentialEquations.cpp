#include <iostream>
#include <clocale>
#include <iomanip>
#include "Solver.h"
#include "DifferentialEquations.h"
using namespace std;

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "��������� ������� ������ ���� ��� ������������� ����������������� ��������� ������� �������\n";
    cout << "������� 13:\n";
    cout << "y' = -y + y^2 + 1\n";

    double h;
    int N;
    double y0;
    double x0 = 0;

    cout << "������� �������� ������� h: ";
    cin >> h;
    cout << "\n";

    cout << "������� ����� ����� ������ �� �0 (N): ";
    cin >> N;
    cout << "\n";

    cout << "������� y(0): ";
    cin >> y0;
    cout << "\n";
    Solver* solver = new Solver(0, y0, h, N);

    for (int k = -2; k <= N; ++k)
    {
        double x = x0 + k * h;
        cout << setprecision(14) << "x: " << x << "\n";

        double preciseSolution = solver->presiseSolution(x);
        cout << setprecision(14) << "������ �������: " << preciseSolution << "\n";

        cout << "����� �������: \n";

        double y = solver->Taylor(x);
        printValue(x, y);
        printPrecision(solver, y, preciseSolution);

        if (k > 2)
        {
            cout << "����� ������: \n";

            double yAdams = solver->Adams(k); 
            printValue(x, yAdams);
            if (k == N) printPrecision(solver, yAdams, preciseSolution);
        }

        if (k > 0)
        {
            cout << "����� �����-�����: \n";

            double yRungeKutta = solver->RungeKutta(k);
            printValue(x, yRungeKutta);
            if (k == N) printPrecision(solver, yRungeKutta, preciseSolution);

            cout << "����� ������: \n";

            double yEuler = solver->Euler(x, k);
            printValue(x, yEuler);
            if (k == N) printPrecision(solver, yEuler, preciseSolution);

            cout << "����� ������ I: \n";

            double yEulerI = solver->EulerI(x, k);
            printValue(x, yEulerI);
            if (k == N) printPrecision(solver, yEulerI, preciseSolution);

            cout << "����� ������ II: \n";

            double yEulerII = solver->EulerII(x, k);
            printValue(x, yEulerII);
            if (k == N) printPrecision(solver, yEulerII, preciseSolution);
        }
        cout << "\n";
    }
}

void printValue(double x, double yAdams)
{
    cout << setprecision(14) << "�������� � ����� " << x << " ����� " << yAdams << "\n";
}

void printPrecision(Solver* solver, double y, double preciseSolution)
{
    double precision = solver->precision(y, preciseSolution);
    cout << setprecision(14) << "���������� �����������: " << precision << "\n\n";
}
