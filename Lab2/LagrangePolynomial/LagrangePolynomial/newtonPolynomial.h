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
		for (int i = 1; i < n; i++)
		{
			double Ai = makeAi(i);
			for (int k = 0; k < i; k++) 
				Ai *= (x - xValues[k]);
			P += Ai;
		}
		return P;
	}

private:
	double makeAi(int i)
	{
		double Ai = 0;
		for (int j = 0; j <= i; j++)
		{
			double w = 1;
			for (int k = 0; k <= i; k++)
			{
				if (k != j)
					w *= (xValues[j] - xValues[k]);
			}
			Ai += yValues[j] / w;
		}
		return Ai;
	}
};