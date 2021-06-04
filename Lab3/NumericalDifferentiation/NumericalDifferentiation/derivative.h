#pragma once
#include "function.h"
#include <vector>
using namespace std;

/// <summary>
/// Подсчёт производных функции методом численного дифференцирования
/// </summary>
class Derivative
{
private:
	double h;

	/// <summary>
	/// Значения функции в точках х
	/// </summary>
	vector<double> xValues;

	/// <summary>
	/// Значения функции в точках у
	/// </summary>
	vector<double> yValues;

	/// <summary>
	/// Значения подсчитаннных значений первых производных в точках
	/// </summary>
	vector<double> firstDerivatives;

	/// <summary>
	/// Погрешности вычислений значений первых производных
	/// </summary>
	vector<double> fstDrvPrecisions;

	/// <summary>
	/// Значения подсчитаннных значений вторых производных в точках
	/// </summary>
	vector<double> secondDerivatives;

	/// <summary>
	/// Погрешности вычислений значений вторых производных
	/// </summary>
	vector<double> sndDrvPrecisions;

	/// <summary>
	/// Первая производная для средних точек в таблице
	/// </summary>
	/// <returns> Приближённое значение первой производной </returns>
	double fstDeriv(double x)
	{
		return (function(x + h) - function(x - h)) / (2 * h);
	}

	/// <summary>
	/// Первая производная для первой точки в таблице
	/// </summary>
	double fstDerivForFirstPoint(double x)
	{
		return (-3 * function(x) + 4 * function(x + h) - function(x + 2 * h)) / (2 * h);
	}

	/// <summary>
	/// Первая производная для последней точки в таблице
	/// </summary>
	double fstDerivForLastPoint(double x)
	{
		return(3 * function(x) - 4 * function(x - h) + function(x - 2 * h)) / (2 * h);
	}

	/// <summary>
	/// Вторая производная для средних точек в таблице
	/// </summary>
	/// <returns> Приближённое значение второй производной </returns>
	double sndDeriv(double x)
	{
		return (function(x + h) - 2 * function(x) + function(x - h)) / (h * h);
	}

public:
	/// <summary>
	/// Конструктор класса
	/// </summary>
	/// <param name="a"> Левая точка отрезка, на котором вычисляем </param>
	/// <param name="h"> Интервал между значениями х </param>
	/// <param name="m"> Число точек </param>
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
	/// Заполнение таблицы значениями первых производных и погрешностей
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
	/// Заполнение таблицы значениями вторых производных и погрешностей
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
	/// Вывод общей таблицы на экран
	/// </summary>
	void printTable()
	{
		calcFstDerivative();
		calcSndDerivative();
		cout << "Таблица значений функции и её производных:\n";

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