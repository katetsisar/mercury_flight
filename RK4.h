#ifndef GAUSS_RK4_H
#define GAUSS_RK4_H

#include <vector>
#include <string>

#include <iostream>
#include <fstream>

using std::vector;
using std::string;

typedef vector<vector<double>> Matrix;



#include "gauss.h"
#include "SubFunctions.h"

// Интегрирование системы ДУ
//
// Входные параметры
//		differentiate     - функция вычисления производных 
//		dimension    - количество параметров
//		steps        - количество шагов
//		x0           - начальное время
//		dx           - шаг по времени
//		y            - вектор начальных условий
//		exit_signal  - сигнал для выхода (0)
//		step_counter - номер шага
//		result       - результирующая матрица
//

inline void differentiate(int dimension,
    double x0,//начальное время
    vector<double> y,
    vector<double>& dy,
    int& exit_signal,//сигнал для выхода (0)
    unsigned int step_number,//номер шага
    Matrix& result);//функция вычисления производных

void RK4(void (*differentiate)(
    int dimension,
    double x0,
    vector<double> y,
    vector<double>& dy,
    int& exit_signal,
    unsigned int step_number, Matrix& result),
    unsigned int dimension, int steps,
    double x0, double dx, vector<double> y,
    int& exit_signal, int& step_counter, Matrix& result);

#endif //GAUSS_RK4_H
