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

// �������������� ������� ��
//
// ������� ���������
//		differentiate     - ������� ���������� ����������� 
//		dimension    - ���������� ����������
//		steps        - ���������� �����
//		x0           - ��������� �����
//		dx           - ��� �� �������
//		y            - ������ ��������� �������
//		exit_signal  - ������ ��� ������ (0)
//		step_counter - ����� ����
//		result       - �������������� �������
//

inline void differentiate(int dimension,
    double x0,//��������� �����
    vector<double> y,
    vector<double>& dy,
    int& exit_signal,//������ ��� ������ (0)
    unsigned int step_number,//����� ����
    Matrix& result);//������� ���������� �����������

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
