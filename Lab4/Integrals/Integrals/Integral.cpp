#include "Integral.h"
#include <cmath>

const double goldenRatio = (1 + sqrt(5)) / 2;
const double accuracy = 0.000001;

double Integral::absPrecision(double J_h)
{
	return abs(J - J_h);
}

double Integral::countMax(double a, double b, int d)
{
	double x1, x2;
	double max = 0.0;
	switch (d + 1)
	{
	case 1:
	{
		max = (*fun).derivative(b);
		break;
	}
	case 2:
	{
		max = (*fun).sndDerivative(b);
		break;
	}
	case 4:
	{
		max = (*fun).fourthDerivative(b);
		break;
	}
	default:
		break;
	}
	return max;
}

double Integral::teoreticalPrecision(double C, int d)
{
	function<double(double)> deriv;
	double M_d_plus_1 = countMax(A, B, d);
	double s = pow(h, d + 1);
	return C * M_d_plus_1 * (B - A) * pow(h, d + 1);
}

void Integral::print(double J_h, double absPres)
{
	printf("����������� �������� ���������: %.14lf \n", J_h);
	printf("����������� �����������: %.14lf \n", absPres);
}

Integral::Integral(Function &fun, double A, double B, int m)
{
	this->A = A;
	this->B = B;
	h = (B - A) / m;
	this->fun = &fun;
	this->J = fun.integral(B) - fun.integral(A);

	for (int k = 1; k < m; ++k)
	{
		double value = A + h * k;
		y += fun.func(value);
	}

	for (int k = 0; k <= m - 1; ++k)
	{
		p += fun.func(A + h * k + h / 2);
	}

	w = fun.func(A) + fun.func(B);
}

void Integral::leftRectangles()
{
	cout << "������� ����� ���������������: \n";

	double J_h = h * ((*fun).func(A) + y);
	double absPres = absPrecision(J_h);

	print(J_h, absPres);
}

void Integral::rightRectangles()
{
	cout << "������� ������ ���������������: \n";

	double J_h = h * (y + (*fun).func(B));
	double absPres = absPrecision(J_h);

	print(J_h, absPres);
}

void Integral::middleRectangles()
{
	cout << "������� ������� ���������������: \n";

	double J_h = h * p;
	double absPres = absPrecision(J_h);

	print(J_h, absPres);
}

void Integral::trapezes()
{
	cout << "������� ��������: \n";

	double J_h = h * (w / 2 + y);
	double absPres = absPrecision(J_h);

	print(J_h, absPres);
}

void Integral::Simpson()
{
	cout << "������� ��������: \n";

	double J_h = h * (p * 2.0 / 3 + w * 1.0 / 6 + y * 1.0 / 3);
	double absPres = absPrecision(J_h);

	print(J_h, absPres);
}

void Integral::countIntegral()
{
	printf("������ �������� ���������: %.14lf \n", J);

	leftRectangles();
	rightRectangles();
	middleRectangles();
	trapezes();
	Simpson();
	printf("\n");
}

void Integral::countIntegralWithTeorPrecision()
{
	printf("������ �������� ���������: %.14lf \n", J);

	leftRectangles();
	double teorPres = teoreticalPrecision(0.5, 0);
	printf("������������� �����������: %.14lf \n", teorPres);

	rightRectangles();
	printf("������������� �����������: %.14lf \n", teorPres);

	middleRectangles();
	teorPres = teoreticalPrecision(1.0 / 24, 1);
	printf("������������� �����������: %.14lf \n", teorPres);

	trapezes();
	teorPres = teoreticalPrecision(1.0 / 12, 1);
	printf("������������� �����������: %.14lf \n", teorPres);

	Simpson();
	teorPres = teoreticalPrecision(1.0 / 2880, 3);
	printf("������������� �����������: %.14lf \n", teorPres);

	ZValues.clear();
}
