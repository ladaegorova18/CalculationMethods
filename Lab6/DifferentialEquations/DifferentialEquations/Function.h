#pragma once
#include <cmath>

class Function
{
private:
	double x0;
	double y0;
	double fstDeriv;
	double sndDeriv;
	double trdDeriv;
	double fourthDeriv;
	double fifthDeriv;
	double sixthDeriv;
	double seventhDeriv;

public:
	Function(double x0, double y0)
	{
		this->x0 = x0;
		this->y0 = y0;

		fstDeriv = -y0 + pow(y0, 2) + 1; //
		sndDeriv = -fstDeriv + 4 * fstDeriv * y0; // 0
		trdDeriv = -sndDeriv + 4 * pow(fstDeriv, 2) + 4 * y0 * sndDeriv; // 
		fourthDeriv = -trdDeriv + 12 * fstDeriv * sndDeriv + 4 * y0 * trdDeriv; // 0
		fifthDeriv = -fourthDeriv + 12 * pow(sndDeriv, 2) + 16 * fstDeriv * trdDeriv + 4 * y0 * fourthDeriv; // 
		sixthDeriv = -fifthDeriv + 20 * fstDeriv * fourthDeriv + 4 * y0 * fifthDeriv; // 0
		seventhDeriv = 4 * y0 * sixthDeriv - sixthDeriv + 40 * pow(trdDeriv, 2) 
			+ 24 * fifthDeriv * fstDeriv + 60 * fourthDeriv * sndDeriv; //

		fstDeriv = 0.875;
		sndDeriv = 0;
		trdDeriv = 3.0625;
		fourthDeriv = 0;
		fifthDeriv = 42.875;
		sixthDeriv = 0;
		seventhDeriv = 1275.53125;
	}

	double function(double x, double y)
	{
		return -y + pow(y, 2) + 1;
	}

	double iDerivative(int i)
	{
		switch (i)
		{
		case 0:
		{
			return y0;
			break;
		}
		case 1:
		{
			return fstDeriv;
			break;
		}
		case 2:
		{
			return sndDeriv;
			break;
		}
		case 3:
		{
			return trdDeriv;
			break;
		}
		case 4:
		{
			return fourthDeriv;
			break;
		}
		case 5:
		{
			return fifthDeriv;
			break;
		}
		case 6:
		{
			return sixthDeriv;
			break;
		}
		case 7:
		{
			return seventhDeriv;
			break;
		}
		default:
			break;
	}
	}
};

