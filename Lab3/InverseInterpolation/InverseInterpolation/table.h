#pragma once
#include <list>
#include <tuple>

std::list<std::tuple<double, double>> makeTable(double a, double b, int m);

std::list<std::tuple<double, double>> sortTable(std::list<std::tuple<double, double>> table, double x);

std::list<std::tuple<double, double>> modifiedSort(std::list<std::tuple<double, double>> table, double F);

std::list<std::tuple<double, double>> swapTable(std::list<std::tuple<double, double>> table);