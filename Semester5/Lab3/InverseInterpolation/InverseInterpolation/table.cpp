#include "table.h"
#include "function.h"
#include <cmath>
#include <iostream>
#include <cctype>
using namespace std;

list<tuple<double, double>> makeTable(double a, double b, int m)
{
	list<tuple<double, double>> table;
	double part = (b - a) / m;
	printf("Zk \t\t\t f(Zk)\n");
	double x = a;
	for (int i = 0; i < m + 1; ++i)
	{
		double y = function(x);
		tuple<double, double> result(x, y);
		table.push_back(result);
		printf("(%.14lf | %.14lf)\n", x, y);
		x += part;
	}
	return table;
}

struct Sort
{
	double X;

public:
	Sort(double x)
	{
		X = x;
	}
	bool operator() (const tuple<double, double> first, const tuple<double, double> second)
	{
		double xFst = 0;
		double xSnd = 0;
		tie(xFst, ignore) = first;
		tie(xSnd, ignore) = second;
		double fstPrecision = abs(xFst - X);
		double sndPrecision = abs(xSnd - X);
		return fstPrecision < sndPrecision;
	}
};

list<tuple<double, double>> sortTable(list<tuple<double, double>> table, double x)
{
	table.sort(Sort(x));
	for (auto couple : table)
	{
		double x = 0;
		double y = 0;
		tie(x, y) = couple;
		printf("(%.14lf | %.14lf)\n", x, y);
	}
	return table;
}

std::list<std::tuple<double, double>> swapTable(std::list<std::tuple<double, double>> table)
{
	for (auto couple = table.begin(); couple != table.end(); ++couple)
	{
		double x = 0;
		double y = 0;
		tie(x, y) = *couple;
		*couple = tuple<double, double>(y, x);
	}
	return table;
}

list<tuple<double, double>> modifiedSort(list<tuple<double, double>> table, double F)
{
	double prevX = 0;
	double currX = 0;
	double c = 0;
	for (auto couple : table)
	{
		double y = 0;
		tie(currX, y) = couple;
		if (F <= y)
			prevX = currX;
		else
		{
			c = (prevX + currX) / 2;
			break;
		}
	}
	return sortTable(table, c);
}