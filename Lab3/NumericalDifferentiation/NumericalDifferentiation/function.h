#pragma once
#include <cmath>

/// <summary>
/// �������, ��� ������� ���� �������� ������ � ������ ����������� � ������
/// </summary>
double function(double x);

/// <summary>
/// ������ ������ ����������� ��� �������� �����������
/// </summary>
double fstDerivative(double x);

/// <summary>
/// ������ ������ ����������� ��� �������� �����������
/// </summary>
double sndDerivative(double x);

/// <summary>
/// ���������� ����������� ��������, ��������� �� ����, � �������� �������
/// </summary>
double precision(double result, double x, double func(double));