#define _USE_MATH_DEFINES

#include "Counter.h"
#include "MelerFunctions.h"
#include <iostream>
#include <limits>
using namespace std;

const double epsilon = 0.0000000000001;

bool is_equal(double x, double y) {
	return std::fabs(x - y) < epsilon;
}

void Counter::countMoments(double A, double B)
{
	moments.push_back(middleRectangles(f->p, A, B));
	moments.push_back(middleRectangles(f->p_1, A, B));
	moments.push_back(middleRectangles(f->p_2, A, B));
	moments.push_back(middleRectangles(f->p_3, A, B));

	for (int i = 0; i < 4; ++i)
	{
		printf("Момент номер %d: %.14lf \n", i, moments[i]);
	}
	printf("\n");
}

double Counter::middleRectangles(double func(double), double A, double B)
{
	double tempM = 100000;
	double h = (B - A) / tempM;
	double p = 0;
	for (int k = 0; k <= tempM; ++k)
	{
		p += func(A + h * k + h / 2);
	}
	return h * p;
}

double Counter::likeGauss()
{
	double denominator = pow(moments[1], 2) - moments[2] * moments[0];
	double a1 = (moments[0] * moments[3] - moments[2] * moments[1]) / denominator;
	double a2 = (pow(moments[2], 2) - moments[3] * moments[1]) / denominator;

	printf("Оргтогональный многочлен w(x): x^2 + %.8lfx + %.8lf\n \n", a1, a2);

	double sqrtD = sqrt(pow(a1, 2) - 4 * a2);
	double x1 = (-a1 + sqrtD) / 2.0;
	double x2 = (-a1 - sqrtD) / 2.0;

	double fraction = 1.0 / (x1 - x2);
	double A1 = fraction * (moments[1] - x2 * moments[0]);
	double A2 = -fraction * (moments[1] - x1 * moments[0]);

	printf("x1: %.14lf \n", x1);
	printf("x2: %.14lf \n", x2);
	printf("A1: %.14lf \n", A1);
	printf("A2: %.14lf \n \n", A2);

	if (x1 != x2)
	{
		printf("Узлы х1 и х2 различны\n");
	}
	if ((A1 > 0 && A2 > 0))
	{
		printf("A1 и A2 положительны\n");
	}
	if (is_equal(A1 + A2, moments[0]))
	{
		printf("Сумма коэффициентов равна нулевому моменту\n");
	}

	return A1 * f->f(x1) + A2 * f->f(x2);
}

Counter::Counter(double A, double B, int m, int N)
{
	this->A = A;
	this->B = B;
	this->m = m;
	this->N = N;
	this->f = new Function();
	this->mehlerF = new MehlerFunctions();
	this->J = f->integral(B) - f->integral(A);
}

void Counter::countLikeGauss()
{
	double integral = 0;
	double tempA = A;
	double part = (B - A) / m;

	countMoments(A, B);
	integral = likeGauss();

	printf("Интеграл по КФ типа Гаусса: %.14lf \n", integral);
	printf("Точное значение интеграла: %.14lf \n", J);

	double precision = f->precision(J, integral);

	printf("Погрешность: %.14lf \n \n", precision);
}

void Counter::countGauss()
{
	double integral = 0;
	double Zk = A;
	double part = (B - A) / m;

	double t1 = -1.0 / sqrt(3);
	double t2 = 1.0 / sqrt(3);

	for (int i = 0; i < m; ++i)
	{
		double x1 = (part / 2) * t1 + Zk + part / 2;
		double x2 = (part / 2) * t2 + Zk + part / 2;

		integral += (part / 2) * (f->func(x1) + f->func(x2));

		Zk += part;
	}

	printf("Интеграл по КФ Гаусса: %.14lf \n", integral);
	printf("Точное значение интеграла: %.14lf \n", J);
	double precision = f->precision(J, integral);

	printf("Погрешность: %.14lf \n \n", precision);
}

void Counter::countMehler()
{
	double integral = 0;
	double mehlerJ = 0;
	double Ai = M_PI / N;
	printf("Точное значение интеграла: %.14lf \n", mehlerJ);

	for (int i = 1; i <= N; ++i)
	{
		double Xk = cos((2.0 * i - 1.0) / (2.0 * N) * M_PI);
		integral += mehlerF->f(Xk);
		printf("X%d: %.14lf \n \n", i, Xk);
		printf("Коэффициент %d: %.14lf \n", i, Ai);
	}
	integral *= Ai;

	printf("Интеграл по КФ Мелера: %.14lf \n", integral);
	printf("Точное значение интеграла: %.14lf \n", mehlerJ);

	double precision = f->precision(mehlerJ, integral);
	printf("Погрешность: %.14lf \n \n", precision);
}
