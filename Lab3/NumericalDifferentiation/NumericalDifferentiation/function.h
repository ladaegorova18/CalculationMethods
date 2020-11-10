#pragma once
#include <cmath>

/// <summary>
/// Функция, для которой ищем значения первой и второй производной в точках
/// </summary>
double function(double x);

/// <summary>
/// Точная первая производная для проверки погрешности
/// </summary>
double fstDerivative(double x);

/// <summary>
/// Точная вторая производная для проверки погрешности
/// </summary>
double sndDerivative(double x);

/// <summary>
/// Вычисление погрешности значения, поданного на вход, и значения функции
/// </summary>
double precision(double result, double x, double func(double));