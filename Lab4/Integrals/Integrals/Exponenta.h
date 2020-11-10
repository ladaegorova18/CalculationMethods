#pragma once
#include "Function.h"

class Exponenta : public Function
{
public:
	Exponenta() {}

	double func(double x) override;

	double integral(double x) override;

	double derivative(double x) override;

	double sndDerivative(double x) override;

	double fourthDerivative(double x) override;
};
