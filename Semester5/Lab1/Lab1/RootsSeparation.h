#pragma once
#ifndef HEADER_H
#define HEADER_H
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <tuple>
#include <iostream>
using namespace std;

list<tuple<double, double>> separation(double A, double B, int N);

#endif // !HEADER_H