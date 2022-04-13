#include <iostream>
#include <locale>
#include "derivative.h"
using namespace std;

/// <summary>
/// Точка входа в программу - ввод данных пользователем и печать результатов
/// </summary>
int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Задача численного дифференцирования:\n";
    cout << "Вариант 4\n";

    int m_plus_1 = 0;
    int m = 0;
    double a = 0;
    double h = 0;

    char key;
    do
    {
        cout << "Введите число отрезков (m + 1):\n";
        cin >> m_plus_1;
        m = m_plus_1 - 1;

        printf("m: %d\n", m);

        cout << "Введите левую границу a:\n";
        cin >> a;

        cout << "Введите интервал h:\n";
        cin >> h;
        while (h <= 0)
        {
            cout << "h должен быть больше 0, попробуйте снова:\n";
            cin >> h;
        }

        auto der = new Derivative(a, h, m);

        der->printTable();

        cout << "Нажмите q, чтобы выйти, или любую клавишу, чтобы продолжить:\n";
        cin >> key;
    } while (key != 'q');
    return 0;
}
