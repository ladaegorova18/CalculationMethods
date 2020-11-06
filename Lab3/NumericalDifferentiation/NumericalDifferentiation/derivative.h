#pragma once
#include "function.h"
#include <vector>
using namespace std;

/// <summary>
/// ������� ����������� ������� ������� ���������� �����������������
/// </summary>
class Derivative
{
private:
	double h;

	/// <summary>
	/// �������� ������� � ������ �
	/// </summary>
	vector<double> xValues;

	/// <summary>
	/// �������� ������� � ������ �
	/// </summary>
	vector<double> yValues;

	/// <summary>
	/// �������� ������������� �������� ������ ����������� � ������
	/// </summary>
	vector<double> firstDerivatives;

	/// <summary>
	/// ����������� ���������� �������� ������ �����������
	/// </summary>
	vector<double> fstDrvPrecisions;

	/// <summary>
	/// �������� ������������� �������� ������ ����������� � ������
	/// </summary>
	vector<double> secondDerivatives;

	/// <summary>
	/// ����������� ���������� �������� ������ �����������
	/// </summary>
	vector<double> sndDrvPrecisions;

	/// <summary>
	/// ������ ����������� ��� ������� ����� � �������
	/// </summary>
	/// <returns> ����������� �������� ������ ����������� </returns>
	double fstDeriv(double x)
	{
		return (function(x + h) - function(x - h)) / (2 * h);
	}

	/// <summary>
	/// ������ ����������� ��� ������ ����� � �������
	/// </summary>
	double fstDerivForFirstPoint(double x)
	{
		return (-3 * function(x) + 4 * function(x + h) - function(x + 2 * h)) / (2 * h);
	}

	/// <summary>
	/// ������ ����������� ��� ��������� ����� � �������
	/// </summary>
	double fstDerivForLastPoint(double x)
	{
		return(3 * function(x) - 4 * function(x - h) + function(x - 2 * h)) / (2 * h);
	}

	/// <summary>
	/// ������ ����������� ��� ������� ����� � �������
	/// </summary>
	/// <returns> ����������� �������� ������ ����������� </returns>
	double sndDeriv(double x)
	{
		return (function(x + h) - 2 * function(x) + function(x - h)) / (h * h);
	}

public:
	/// <summary>
	/// ����������� ������
	/// </summary>
	/// <param name="a"> ����� ����� �������, �� ������� ��������� </param>
	/// <param name="h"> �������� ����� ���������� � </param>
	/// <param name="m"> ����� ����� </param>
	Derivative(double a, double h, int m)
	{
		this->h = h;
		for (double i = 0; i <= m; ++i)
		{
			double value = a + i * h;
			xValues.push_back(value);
			yValues.push_back(function(value));
		}
	}

	/// <summary>
	/// ���������� ������� ���������� ������ ����������� � ������������
	/// </summary>
	void calcFstDerivative()
	{
		firstDerivatives.push_back(fstDerivForFirstPoint(xValues[0]));
		for (int i = 1; i < xValues.size() - 1; ++i)
		{
			firstDerivatives.push_back(fstDeriv(xValues[i]));
		}
		firstDerivatives.push_back(fstDerivForLastPoint(xValues[xValues.size() - 1]));

		for (int i = 0; i < xValues.size(); ++i)
		{
			fstDrvPrecisions.push_back(precision(firstDerivatives[i], xValues[i], fstDerivative));
		}
	}

	/// <summary>
	/// ���������� ������� ���������� ������ ����������� � ������������
	/// </summary>
	void calcSndDerivative()
	{
		secondDerivatives.push_back(0);
		sndDrvPrecisions.push_back(0);
		for (int i = 1; i < xValues.size() - 1; ++i)
		{
			secondDerivatives.push_back(sndDeriv(xValues[i]));
		}

		for (int i = 1; i < xValues.size() - 1; ++i)
		{
			sndDrvPrecisions.push_back(precision(secondDerivatives[i], xValues[i], sndDerivative));
		}
		secondDerivatives.push_back(0);
		sndDrvPrecisions.push_back(0);
	}

	/// <summary>
	/// ����� ����� ������� �� �����
	/// </summary>
	void printTable()
	{
		calcFstDerivative();
		calcSndDerivative();
		cout << "������� �������� ������� � � �����������:\n";

		cout << " \t x \t \t f(x) \t \t f'(x) calc.   |f'(x)-f'(x) calc.| f\"(x) calc. |f\"(x) - f\" calc.|\n";
		for (int i = 0; i < xValues.size(); ++i)
		{
			if (i != 0 || i != xValues.size() - 1)
			{
				printf("%.14lf | %.14lf | %.14lf | %.14lf | %.14lf | %.14lf\n", xValues[i], yValues[i],
					firstDerivatives[i], fstDrvPrecisions[i], secondDerivatives[i], sndDrvPrecisions[i]);
			}
			else
			{
				printf("%.14lf | %.14lf | %.14lf | %.14lf | ---- | ---- |\n", xValues[i], yValues[i],
					firstDerivatives[i], fstDrvPrecisions[i]);
			}
		}
	}
};