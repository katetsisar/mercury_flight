#ifndef GAUSS_ITERSYSTEMSOLVER_H
#define GAUSS_ITERSYSTEMSOLVER_H

#include <iostream>
#include <vector>
#include <string>

//#include "SubFunctions.h"

using std::vector;

typedef vector<vector<double>> Matrix;

// Решение СЛАУ
//
// Входные параметры
//		LU - треугольные матрицы
//		B  - вектор свободных членов
//		T  - массив перестановок строк матрицы A
//		X  - вектор решения системы
//
void solveSystem(
	const Matrix& LU, vector<double>& B, 
	const vector<int>& T, vector<double>& X);

#endif //GAUSS_ITERSYSTEMSOLVER_H
