#pragma once
#include <list>
#include <tuple>

std::list<std::tuple<double, double>> makeTable(double a, double b, int m);

std::list<std::tuple<double, double>> sortTable(std::list<std::tuple<double, double>> table, double x);