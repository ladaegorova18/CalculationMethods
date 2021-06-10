#pragma once
#include <iostream>
#include <clocale>
#include <vector>
using namespace std;

vector<vector<double>> countX(vector<vector<double>> M, vector<vector<double>> b);

vector<vector<double>> countDelta(vector<vector<double>> A, vector<vector<double>> A_);

double norm(vector<vector<double>> A);
