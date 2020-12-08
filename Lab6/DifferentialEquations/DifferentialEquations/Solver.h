#pragma once
#include "Function.h"
#include <vector>
using namespace std;

class Solver
{
private:
	double x0;
	double y0;
	double c;
	double h;
	double adams[5];
	int N;

	vector<double> xValues;

	vector<double> RKvalues;
	vector<double> X;
	vector<double> Y;
	vector<double> EulerValues;
	vector<double> EulerIValues;
	vector<double> EulerIIValues;

	vector<vector<double>> table;

	
	Function* func;

	int factorial(int i)
	{
		return (i == 0 || i == 1) ? 1 : factorial(i - 1) * i;
	}

public:
	Solver(double x0, double y0, double h, int N)
	{
		this->x0 = x0;
		this->y0 = y0;
		this->h = h;
		this->N = N;

		for (int i = 0; i < 5; ++i)
		{
			vector<double> v(N + 2);
			table.push_back(v);
		}

		for (int i = -2; i < N; ++i)
		{
			X.push_back(x0 + h * i);
		}

		for (int i = 0; i <= N; ++i)
		{
			xValues.push_back(x0 + i * h);
		}
		c = (atan((4 * y0 - 1) / sqrt(7)) * 2 - sqrt(7) * x0) / (sqrt(7));

		func = new Function(x0, y0);
		RKvalues.push_back(y0);
		EulerValues.push_back(y0);
		EulerIValues.push_back(y0);
		EulerIIValues.push_back(y0);

	}

	double presiseSolution(double x)
	{
		return 0.25 * (sqrt(7) * tan(0.5 * sqrt(7) * (c + x)) + 1);
	}

	double Taylor(double x)
	{
		double y = 0;
		for (int i = 0; i <= 7; ++i)
		{
			y += func->iDerivative(i) * pow(x - x0, i) / factorial(i);
		}
		if (Y.size() < 5)
			Y.push_back(y);
		return y;
	}

	double Adams(int n)
	{
		n += 1;
		
		for (int i = n - 4; i <= n; ++i)
		{
			table[0][i] = h * func->function(X[i], Y[i]);
		}
		/*for (int i = n - 4; i <= n - 1; ++i)
		{
			table[1][i] = table[0][i + 1] - table[0][i];
		}
		for (int i = n - 4; i <= n - 2; ++i)
		{
			table[2][i] = table[1][i + 1] - table[1][i];
		}
		for (int i = n - 4; i <= n - 3; ++i)
		{
			table[3][i] = table[2][i + 1] - table[2][i];
		}
		table[4][n - 4] = table[3][n - 3] - table[3][n - 4];*/

		for (int j = 4; j > 0; j--)
		{
			for (int k = 0; k < j; k++)
			{
				table[k][5 - j] = table[k + 1][4 - j] - table[k][4 - j];
			}
		}
		double deltayn = table[0][n] + 0.5 * table[1][n - 1] + (5.0 / 12) * table[2][n - 2] + (3.0 / 8) * table[3][n - 3]
			+ (251.0 / 720) * table[4][n - 4];
		double y = Y[n - 1] + deltayn;
		Y.push_back(y);
		return y;
	}

	double RungeKutta(int i)
	{
		double x = xValues[i - 1];
		double y = RKvalues[i - 1];
		double k1 = h * func->function(x, y);
		double k2 = h * func->function(x + h / 2, y + k1 / 2);
		double k3 = h * func->function(x + h / 2, y + k2 / 2);
		double k4 = h * func->function(x + h, y + k3);

		double result = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		RKvalues.push_back(result);
		return result;
	}

	double Euler(double x, int k)
	{
		double y = EulerValues[k - 1] + h * func->function(xValues[k - 1], EulerValues[k - 1]);
		EulerValues.push_back(y);
		return y;
	}

	double EulerI(double x, int k)
	{
		double y_k_plus_half = EulerIValues[k - 1] + (h / 2) * func->function(xValues[k - 1], EulerIValues[k - 1]);
		double y = EulerIValues[k - 1] + h * func->function(xValues[k - 1] + h / 2, y_k_plus_half);
		EulerIValues.push_back(y);
		return y;
	}

	double EulerII(double x, int k)
	{
		double f = func->function(xValues[k - 1], EulerIIValues[k - 1]);
		double Y_k = EulerIIValues[k - 1] + h * f;
		double y = EulerIIValues[k - 1] + (h / 2) * (f + func->function(xValues[k], Y_k));
		EulerIIValues.push_back(y);
		return y;
	}

	double precision(double result, double exactValue)
	{
		return abs(result - exactValue);
	}
};



//for (int i = n - 4; i <= n; ++i)
//{
//	table[0][i] = h * func->equation(X[i], Y[i]);
//}
//for (int i = n - 4; i <= n - 1; ++i)
//{
//	table[1][i] = table[0][i + 1] - table[0][i];
//}
//for (int i = n - 4; i <= n - 2; ++i)
//{
//	table[2][i] = table[1][i + 1] - table[1][i];
//}
//for (int i = n - 4; i <= n - 3; ++i)
//{
//	table[3][i] = table[2][i + 1] - table[2][i];
//}
//table[4][n - 4] = table[3][n - 3] - table[3][n - 4];



/*for (int i = n - 4; i <= n - 1; ++i)
{
	table[1][i] = table[0][i + 1] - table[0][i];
}
for (int i = n - 4; i <= n - 2; ++i)
{
	table[2][i] = table[1][i + 1] - table[1][i];
}
for (int i = n - 4; i <= n - 3; ++i)
{
	table[3][i] = table[2][i + 1] - table[2][i];
}
table[4].push_back(table[3][n - 3] - table[3][n - 4]);*/

//double deltayn = table[0][n] + 0.5 * table[1][n - 1] + (5.0 / 12) * table[2][n - 2] + (3.0 / 8) * table[3][n - 3]
//	+ (251.0 / 720) * table[4][n - 4];
//double y = Y[n - 1] + deltayn;
//Y.push_back(y);
//return y;



//double delta = 0;
//vector<double> result(N + 1);
//vector<vector<double>> table;
//
//for (int i = 0; i < N + 3; ++i)
//{
//	vector<double> v(5);
//	table.push_back(v);
//}
//
//for (int i = 0; i < 5; ++i)
//{
//	table[0][i] = h * func->function(X[i], Y[i]);
//}
//
//for (int j = 4; j > 0; j--)
//{
//	for (int k = 0; k < j; k++)
//	{
//		table[k][5 - j] = table[k + 1][4 - j] - table[k][4 - j];
//	}
//}
//for (int k = 3; k < N + 1; k++)
//{
//	delta = table[k + 1][0] + (table[k][1]) / 2.0 + (table[k - 1][2] * 5) / 12.0 + (table[k - 2][3] * 3) / 8.0 + (table[k - 3][4] * 251) / 720.0;
//
//	result[k] = result[k - 1] + delta;
//
//	table[k + 2][0] = h * func->function(X[k], result[k]);
//	table[k + 1][1] = table[k + 2][0] - table[k + 1][0];
//	table[k][2] = table[k + 1][1] - table[k][1];
//	table[k - 1][3] = table[k][2] - table[k - 1][2];
//	table[k - 2][4] = table[k - 1][3] - table[k - 2][3];
//}
//return result;