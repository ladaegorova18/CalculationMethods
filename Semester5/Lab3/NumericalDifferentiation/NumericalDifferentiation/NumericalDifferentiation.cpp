#include <iostream>
#include <locale>
#include "derivative.h"
using namespace std;

/// <summary>
/// ����� ����� � ��������� - ���� ������ ������������� � ������ �����������
/// </summary>
int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "������ ���������� �����������������:\n";
    cout << "������� 4\n";

    int m_plus_1 = 0;
    int m = 0;
    double a = 0;
    double h = 0;

    char key;
    do
    {
        cout << "������� ����� �������� (m + 1):\n";
        cin >> m_plus_1;
        m = m_plus_1 - 1;

        printf("m: %d\n", m);

        cout << "������� ����� ������� a:\n";
        cin >> a;

        cout << "������� �������� h:\n";
        cin >> h;
        while (h <= 0)
        {
            cout << "h ������ ���� ������ 0, ���������� �����:\n";
            cin >> h;
        }

        auto der = new Derivative(a, h, m);

        der->printTable();

        cout << "������� q, ����� �����, ��� ����� �������, ����� ����������:\n";
        cin >> key;
    } while (key != 'q');
    return 0;
}
