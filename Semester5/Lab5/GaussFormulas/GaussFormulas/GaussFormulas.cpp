#include <iostream>
#include <clocale>
#include <vector>
#include "Counter.h"
using namespace std;

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "����������� ���������� ���������� ��� ������ �� ����\n";
    cout << "������� 13:\n";

    char key = 'w';
    do
    {
        double A;
        double B;
        int m;
        int N = 2;
        cout << "������� A:\n";
        cin >> A;           // A = 0
        cout << "������� B:\n";
        cin >> B;           // B = 1 
        cout << "������� ����� ����������� ������� m:\n";
        cin >> m;           // m = 100
        cout << "������� ����� ����� N ��� �� ������:\n";
        cin >> N;           // N = 6

        Counter* counter = new Counter(A, B, m, N);
        counter->countGauss();
        counter->countLikeGauss();
        counter->countMehler();

        cout << "������� q, ����� �����, ����� ������ �������, ����� ����������:\n";
        cin >> key;
    } while (key != 'q');
}
