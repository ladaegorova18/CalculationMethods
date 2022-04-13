#pragma once
#include <vector>
#include <iostream>
#include <functional>
#include "Function.h"
#include "PolynomialFirstDgr.h"
#include "PolynomialThirdDegree.h"
using namespace std;

class Integral
{
private:
	
	Function *fun;
	double A;
	double B;
	double h;
	double J;
	double y = 0;
	double p = 0;
	double w;
	double m = 0;

	vector<double> ZValues;

	double absPrecision(double J_h);

	double teoreticalPrecision(double C, int d);

	double countMax(double a, double b, int d);

	void print(double J_h, double absPres);

public:
	Integral(Function &fun, double A, double B, int m);

	void leftRectangles();

	void rightRectangles();

	void middleRectangles();

	void trapezes();

	void Simpson();

	void countIntegral();

	void countIntegralWithTeorPrecision();
};

