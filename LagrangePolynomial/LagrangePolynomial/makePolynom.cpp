#include "makePolynom.h"
#include <vector>
#include <tuple>

double polynomial(list<tuple<double, double>> table, int n, double x)
{
	vector<double> xValues;
	vector<double> yValues;

	auto it = table.begin();
	for (int k = 0; k < n + 1; ++k)
	{
		double x = 0;
		double y = 0;
		tie(x, y) = *it;
		++it;
		xValues.push_back(x);
		yValues.push_back(y);
	}

	double P = 0.0;

	for (int k = 0; k < n + 1; ++k)
	{
		double  lk = 1.0;
		for (int j = 0; j < n + 1; ++j)
		{
			if (k != j)
			{
				lk *= (x - xValues[j]) / (xValues[k] - xValues[j]);
			}
		}
		P = P + yValues[k] * lk;
	}
	return P;
}