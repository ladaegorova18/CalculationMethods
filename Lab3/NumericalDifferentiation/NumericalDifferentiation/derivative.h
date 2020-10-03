#pragma once
#include "function.h"
#include <vector>
using namespace std;

class Derivative
{
private:
	double h;

	vector<double> xValues;
	vector<double> yValues;

	vector<double> firstDerivatives;
	vector<double> fstDrvPrecisions;

	vector<double> secondDerivatives;
	vector<double> sndDrvPrecisions;

	double fstDeriv(double x)
	{
		return (function(x + h) - function(x - h)) / (2 * h);
	}

	double fstDerivForFirstPoint(double x)
	{
		return (-3 * function(x) + 4 * function(x + h) - function(x + 2 * h)) / (2 * h);
	}

	double fstDerivForLastPoint(double x)
	{
		return(3 * function(x) - 4 * function(x - h) + function(x - 2 * h)) / (2 * h);
	}

	double sndDeriv(double x)
	{
		return (function(x + h) - 2 * function(x) + function(x - h)) / (h * h);
	}

public:
	Derivative(double a, double h, int m)
	{
		this->h = h;
		for (double i = 0; i < m; ++i)
		{
			double value = a + i * h;
			xValues.push_back(value);
			yValues.push_back(function(value));
		}
	}

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

	void printTable()
	{
		calcFstDerivative();
		calcSndDerivative();
		cout << "Таблица значений функции и её производных:\n";

		cout << " \t x \t \t f(x) \t \t f'(x) calc.\t|f'(x)-f'(x) calc.|\t f\"(x) calc. \t |f\"(x) - f\" calc.|\n";
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