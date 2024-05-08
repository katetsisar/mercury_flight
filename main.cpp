#include "RK4.h"
#include "Modyfied_Newton.h"
#include "SubFunctions.h"


//���������� ���������
//declaring global variables
//�������� ��� ������
//extern double A0 = 0.00573449023;
extern double A0 = 0.0059361;
//extern double V0 = 29410.51510721725;//��������� ���� ��������, ����� ����� �������� �����, �/c
extern double V0 = 29800;
extern double U0 = 0;//��������� ��� ��������,  �/c
//extern double R0 = 152098232000;//������ ������ ������, �
extern double R0 = 149600000000;
//extern double VM= 47000.36;//����������� �������� ��������,  �/c
extern double VM = 38851;
//extern double RM= 69817445000;//������ ������ ��������, �
extern double RM = 69817400000;
extern double Totn = 8.3e-4;//��������� �� ���� ����������
extern double FiZ = 0.0;//��������� �������� ���� �����
extern double Qotn = 1.29e-3;
extern double tabs = 0;
extern double steps = 1000;
//��������� ����������
extern double t1 = 200;
extern double PU = -621;
extern double PV = -911;
extern double PR = -0.00001;
//extern double PR = 0;


vector<double> split(string str, string delimiter)
{
	vector<string> v;
	if (!str.empty()) {
		int start = 0;
		do {
			// Find the index of occurrence
			int idx = str.find(delimiter, start);
			if (idx == string::npos) {
				break;
			}

			// If found add the substring till that
			// occurrence in the vector
			int length = idx - start;
			v.push_back(str.substr(start, length));
			start += (length + delimiter.size());
		} while (true);
		v.push_back(str.substr(start));
	}
	vector<double> vd(v.size()-1);
	for (int i = 0;i < v.size()-1;i++) {
		vd[i] = std::stod(v[i]);
	}
	return vd;
}

//������� ���� ������ ��� ��4 � ���������
int main() {
	int dimension = 4;
	//����� �������� ���� pr=0(
	vector<double> x = { PU / (PU),PV / (PV),PR / (PR),t1 / (t1) };//���������� ��������� ����������
	vector<double> residuals(4);
	int iter = 0;
    newtonModyfied(&getResidual, dimension, x, residuals, iter);
	std::cout << "\nafter calculations shooting params x=" << x[0] << ";" << x[1] << ";" << x[2] << ";" << x[3] << "\n";
	std::cout << "\nfinal answer: " << x[0] * PU << ';' << x[1] * PV << ';' << x[2] * PR << ';' << x[3] * 200;
	/*
	//result graph rk4 written in "rk4_new.csv"
	//read file result
	std::ifstream file;
	file.open("rk4_new.csv");
	string line;
	std::getline(file, line);
	vector<double> r(13);
	Matrix mBuild(1001, vector<double>(13));
	int j = 0;
	while (std::getline(file, line)) {
		//std::cout << line << std::endl;
		r = split(line, ";");//������ �������� �� ������
		mBuild[j] = r;
		j++;
	}
	int n = 1000;
	std::vector<double> tb(n), ub(n);
	for (int i = 0;i < n;i++) {
		tb[i] = mBuild[i][0];
		ub[i] = mBuild[i][1];
	}
	std::cout << "tb = ";
	for (int i = 0;i < n;i++) {
		std::cout << tb[i] << ";";
	}
	std::cout << std::endl;
	std::cout << "ub = ";
	for (int i = 0;i < n;i++) {
		std::cout << ub[i] << ";";
	}
	std::cout << std::endl;
	*/
	return 0;
}