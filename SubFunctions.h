#ifndef SUBFUNCTIONS_H
#define SUBFUNCTIONS_H

#include <vector>
#include <cmath>
#include "gauss.h"
#include "RK4.h"

using std::vector;
typedef vector<vector<double>> Matrix;

extern double A0;
extern double V0;//��������� ���� ��������, ����� ����� �������� �����, �/c
extern double U0;//��������� ��� ��������,  �/c
extern double R0;//������ ������ ������, �
extern double Totn;//��������� �� ���� ����������(�� �����)
extern double FiZ;//��������� �������� ���� �����
extern double Qotn;//���� ���� ��� � ������� ���������
extern double tabs;
extern double steps;
extern double VM;//����������� �������� ��������,  �/c
extern double RM;//������ ������ ��������, �
extern double t1;
extern double PU;
extern double PV;
extern double PR;


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

inline void differentiate(int dimension,
    double x0,//��������� �����
    vector<double> y,
    vector<double>& dy,
    int& exit_signal,//������ ��� ������ (0)
    unsigned int step_number,//����� ����
    Matrix& result);//������� ���������� �����������

void getResidual(
    const vector<double>& x,
    vector<double>& residuals);


//
double getEuclideanNorm(vector<double>& residuals);

double getLocalNorm(vector<double>&  residuals, Matrix& A);


#endif //GAUSS_GAUSS_H