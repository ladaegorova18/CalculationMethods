#pragma once
#include "Function.h"

class PolynomialFirstDgr : public Function
{
public:
	PolynomialFirstDgr() {}

	double func(double x) override;

	double integral(double x) override;

	double derivative(double x) override;

	double sndDerivative(double x) override;

	double fourthDerivative(double x) override;
};

