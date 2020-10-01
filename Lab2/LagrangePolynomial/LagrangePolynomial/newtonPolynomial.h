#pragma once
#include <list>
#include <tuple>
#include <vector>
using namespace std;

class Newton
{
	vector<double> xValues;
	vector<double> yValues;

	int n = 0;
	double step = 0;

public:
	Newton(list<tuple<double, double>> table, int n)
	{
		auto it = table.begin();
		for (int k = 0; k < n; ++k)
		{
			double x = 0;
			double y = 0;
			tie(x, y) = *it;
			++it;
			xValues.push_back(x);
			yValues.push_back(y);
		}
		this->n = n;
		step = xValues[1] - xValues[0];
	}

	double count(double x)
	{
		double P = yValues[0];

		for (int j = 1; j < n; j++)
		{
			double something = makeAk(j);
			for (int k = 0; k < j; k++)
				something *= x - xValues[k];
			P += something;
		}
		return P;
	}

private:
	unsigned int factorial(int i)
	{
		if (i == 0)
			return 1;
		return i * factorial(i - 1);
	}

	double makeAk(int j)
	{
		if (j == 0)
			return yValues[0];
		double delta = Delta(j, j);
		return delta / (factorial(j) * pow(step, j));
	}

	double Delta(int p, int i)
	{
		if (p == 1)
			return yValues[i] - yValues[i - 1];
		return Delta(p - 1, i) - Delta(p - 1, i - 1);
	}
};