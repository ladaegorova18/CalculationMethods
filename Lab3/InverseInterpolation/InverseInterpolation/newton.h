//#pragma once
//#include <list>
//#include <tuple>
//#include <vector>
//using namespace std;
//
//class Newton
//{
//	vector<double> xValues;
//	vector<double> yValues;
//
//	int n = 0;
//	double step = 0;
//
//public:
//	Newton(list<tuple<double, double>> table, int n)
//	{
//		auto it = table.begin();
//		for (int k = 0; k < n; ++k)
//		{
//			double x = 0;
//			double y = 0;
//			tie(x, y) = *it;
//			++it;
//			xValues.push_back(x);
//			yValues.push_back(y);
//		}
//		this->n = n;
//		step = xValues[1] - xValues[0];
//	}
//
//	double count(double x)
//	{
//		double P = yValues[0];
//
//		for (int j = 1; j < n; j++)
//		{
//			double something = makeAk(j);
//			for (int k = 0; k < j; k++)
//				something *= x - xValues[k];
//			P += something;
//		}
//		return P;
//	}
//
//	double bisection(double a, double b, double epsilon)
//	{
//		int N = 0;
//		double X0 = (a + b) / 2;
//
//		do
//		{
//			double c = (a + b) / 2;
//			if (count(a) * count(c) <= 0)
//			{
//				b = c;
//			}
//			else
//			{
//				a = c;
//			}
//			++N;
//		} while (b - a > 2 * epsilon);
//
//		double X = (a + b) / 2;
//		return X;
//	}
//
//private:
//	unsigned int factorial(int i)
//	{
//		if (i == 0)
//			return 1;
//		return i * factorial(i - 1);
//	}
//
//	double makeAk(int j)
//	{
//		if (j == 0)
//			return yValues[0];
//		double delta = Delta(j, j);
//		return delta / (factorial(j) * pow(step, j));
//	}
//
//	double Delta(int p, int i)
//	{
//		if (p == 1)
//			return yValues[i] - yValues[i - 1];
//		return Delta(p - 1, i) - Delta(p - 1, i - 1);
//	}
//};

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

	double bisection(double a, double b, double epsilon, double F)
	{
		int N = 0;
		double X0 = (a + b) / 2;

		do
		{
			double c = (a + b) / 2;
			if ((count(a) - F) * (count(c) - F) <= 0)
			{
				b = c;
			}
			else
			{
				a = c;
			}
			++N;
		} while (b - a > 2 * epsilon);

		double X = (a + b) / 2;
		return X;
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