#pragma once
#include "Function.h"
#include <vector>
#include <iostream>
#include <iomanip>
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

		for (int i = -2; i <= N; ++i)
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
		return 0.25 * (sqrt(7) * tan(0.5 * sqrt(7) * x) + 1);
	}

	double Taylor(double x)
	{
		double y = 0;
		for (int i = 0; i <= 7; ++i)
		{
			y += func->iDerivative(i) * pow(x - x0, i) / factorial(i);
		}
		Y.push_back(y);
		return y;
		//return 0.25 + 7 * x / 8 + 49 * pow(x, 3) / 96 + 343 * pow(x, 5) / 960 + 5831 * pow(x, 7) / 23040;
	}

	vector<double> Adams(int n)
	{
		vector<double> q, dq, dSqrq, dCubeq, d4q, results;
		results.push_back(y0);

		for (int i = 0; i < 5; i++) {
			q.push_back(h * func->function(X[i], Y[i]));
		}
		for (int i = 0; i < 4; i++) {
			dq.push_back(q[i + 1] - q[i]);
		}
		for (int i = 0; i < 3; i++) {
			dSqrq.push_back(dq[i + 1] - dq[i]);
		}
		for (int i = 0; i < 2; i++) {
			dCubeq.push_back(dSqrq[i + 1] - dSqrq[i]);
		}
		for (int i = 0; i < 1; i++) {
			d4q.push_back(dCubeq[i + 1] - dCubeq[i]);
		}

		for (int i = 5; i <= n + 2; i++) 
		{
			Y[i] = (Y[i - 1] + q[i - 1] + dq[i - 2] / 2.0 + 5 * dSqrq[i - 3] / 12.0 
				+ 3 * dCubeq[i - 4] / 8.0 + 251 * d4q[i - 5] / 720.0);
			q.push_back(h * func->function(X[i], Y[i]));
			dq.push_back(q[i] - q[i - 1]);
			dSqrq.push_back(dq[i - 1] - dq[i - 2]);
			dCubeq.push_back(dSqrq[i - 2] - dSqrq[i - 3]);
			d4q.push_back(dCubeq[i - 3] - dCubeq[i - 4]);
		}

		return Y;
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