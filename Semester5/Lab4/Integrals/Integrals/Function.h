#pragma once

class Function
{
public:
	virtual double func(double x)
	{
		return 0.0;
	}

	virtual double integral(double x)
	{
		return 0.0;
	}

	virtual double derivative(double x)
	{
		return 0.0;
	}

	virtual double sndDerivative(double x)
	{
		return 0.0;
	}

	virtual double fourthDerivative(double x)
	{
		return 0.0;
	}
};

