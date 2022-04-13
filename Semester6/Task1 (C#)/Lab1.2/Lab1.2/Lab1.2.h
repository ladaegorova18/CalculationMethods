#pragma once
#include <iostream>
#include <clocale>
#include <vector>
#include <iomanip>
using namespace std;

void printMatrix(vector<vector<double>> A);

void printResult(vector<double> X, string name);

vector<double> GaussMethod(vector<vector<double>> A);

vector<double> gaussWithMainElement(vector<vector<double>> A);

vector<vector<double>> ForvardMove(vector<vector<double>> A, bool chooseMainElement);

vector<double> ReverseMove(vector<vector<double>> A);

double LUdecomposition(size_t n, vector<vector<double>> A);

vector<vector<double>> Jordan(vector<vector<double>> A);

vector<vector<double>> multiply(vector<vector<double>> revA, vector<vector<double>> A);

vector<vector<double>> reverseMatrix(vector<vector<double>> A, size_t n);

vector<vector<double>> chooseMain(vector<vector<double>> A, int k);

void printPrecision(vector<vector<double>> A, vector<double> X);