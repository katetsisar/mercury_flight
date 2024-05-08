#ifndef GAUSS_MODYFIED_NEWTON_H
#define GAUSS_MODYFIED_NEWTON_H

#include <algorithm>
#include <vector>
#include <string>

#include <iostream>
#include <fstream>

#include "SubFunctions.h"
#include "gauss.h"
#include "iterSystemSolver.h"

using std::vector;
using std::string;


// Решение системы нелинейных алгебраических уравнений
//
// Входные параметры
//		getResidual - функция нахождения невязок
//		dimension   - количество параметров пристрелки
//		x           - параметры пристрелки
//		residuals   - невязки (считаются внутри)
//		iteration   - номер итерации
//


void LUP(const Matrix& A, Matrix& LU, vector<int>& T, string& error);

void newtonModyfied(
    void (*getResidual) (
        const vector<double>& ,
        vector<double>& ),
    unsigned int dimension,
    vector<double>& x,
    vector<double>& residuals,
    int& iteration);

#endif //GAUSS_CLASSIC_NEWTON_H
