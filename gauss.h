#ifndef GAUSS_GAUSS_H
#define GAUSS_GAUSS_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using std::vector;
using std::string;

typedef vector<vector<double>> Matrix;

// ���������� ������� A � ������������ ���� ����������� ������ 
//
// ������� ���������
//		A     - ������� ������� (������� �����)
//		LU    - ����������� �������
//		T     - ������ ������������ ����� ������� A
//              (��������� ��������� ���������������� ������� �� ��������) 
//		error - ��������� ������
//
void LUP(const Matrix &A, Matrix &LU, vector<int> &T, string& error);

#endif //GAUSS_GAUSS_H