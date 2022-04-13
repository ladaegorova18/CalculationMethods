#pragma once
#include <vector>
#include "Function.h"
#include "MelerFunctions.h"
using namespace std;

class Counter
{
private:
	double A;
	double B;
	double J;
	int m;
	int N;
	Function* f;
	MehlerFunctions* mehlerF;

	vector<double> moments;
	void countMoments(double A, double B);
	double middleRectangles(double func(double), double A, double B);
	double likeGauss();

public:
	Counter(double A, double B, int m, int N);
	void countGauss();
	void countLikeGauss();
	void countMehler();
};

