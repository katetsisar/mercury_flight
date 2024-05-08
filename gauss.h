#ifndef GAUSS_GAUSS_H
#define GAUSS_GAUSS_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using std::vector;
using std::string;

typedef vector<vector<double>> Matrix;

// –азложение матрицы A в произведение двух треугольных матриц 
//
// ¬ходные параметры
//		A     - исходн€ матрица (матрица якоби)
//		LU    - треугольные матрицы
//		T     - массив перестановок строк матрицы A
//              (единична€ диагональ нижнетреугольной матрицы не хранитс€) 
//		error - индикатор ошибки
//
void LUP(const Matrix &A, Matrix &LU, vector<int> &T, string& error);

#endif //GAUSS_GAUSS_H