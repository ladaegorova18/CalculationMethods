#include <iostream>
#include <clocale>
#include <vector>
#include "Counter.h"
using namespace std;

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Приближённое вычисление интегралов при помощи КФ НАСТ\n";
    cout << "Вариант 13:\n";

    char key = 'w';
    do
    {
        double A;
        double B;
        int m;
        int N = 2;
        cout << "Введите A:\n";
        cin >> A;           // A = 0
        cout << "Введите B:\n";
        cin >> B;           // B = 1 
        cout << "Введите число промежутков деления m:\n";
        cin >> m;           // m = 100
        cout << "Введите число узлов N для КФ Мелера:\n";
        cin >> N;           // N = 6

        Counter* counter = new Counter(A, B, m, N);
        counter->countGauss();
        counter->countLikeGauss();
        counter->countMehler();

        cout << "Нажмите q, чтобы выйти, любую другую клавишу, чтобы продолжить:\n";
        cin >> key;
    } while (key != 'q');
}
