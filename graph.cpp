#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using std::vector;
using std::string;
typedef vector<vector<double>> Matrix;

#include "graph.h"

//������ ������� �� ������ ��� ������� � ������� ������ �������
//U(t); V(t); R(t), Fi(t); Teta(t), PsiU(t), PsiV(t), PsiR(t), y(x)
//����� ��������� 1000����� �� ���������� ����������� �������

namespace plt = matplotlibcpp;

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

int main()
{
    //result graph rk4 written in "rk4_new.csv"
	//read file result
	std::ifstream file;
	file.open("rk4_new.csv");
	string line;
	std::getline(file, line);
	vector<double> r(13);
	Matrix mb(1001, vector<double>(13));
	int j = 0;
    while (std::getline(file, line)) {
        //std::cout << line << std::endl;
        r = split(line, ";");//вектор значений из строки
        mb[j] = r;
        j++;
    }
    int n = 1001;
    std::vector<double> tb(n),ub(n), vb(n), rb(n), fib(n), pub(n), pvb(n), prb(n), pfib(n), hpb(n), sintb(n), costb(n), tetab(n),xb(n),yb(n);
    //std::cout << mb[0].size()<<"\n";
    
    std::cout<<mb[0][1] << "\n";
    
    /*
    std::cout << mb[0][1] << "\n";
    std::cout << mb[0][2] << "\n";
    std::cout << mb[0][3] << "\n";
    std::cout << mb[0][4] << "\n";
    std::cout << mb[0][5] << "\n";
    std::cout << mb[0][6] << "\n";
    std::cout << mb[0][7] << "\n";
    std::cout << mb[0][8] << "\n";
    std::cout << mb[0][9] << "\n";
    std::cout << mb[0][10] << "\n";
    std::cout << mb[0][11] << "\n";
    std::cout << mb[0][12] << "\n";
    */

    //std::cout << mb.size() << "\n";
    //НУ заплдняем не из матрицы

    ub[0] = 0;
    vb[0] = 29800;
    rb[0] = 149600000000;
    fib[0] = 0;
    pub[0] = -411.35;
    pvb[0] = -361.438;
    prb[0] = -0.000155226;
    pfib[0] = 0;
    hpb[0] = -0.545506;
    sintb[0] = -0.751212;
    costb[0] = -0.660061;
    tb[0] = 0;
    xb[0] = rb[0] * cos(fib[0]);
    yb[0]= rb[0] * sin(fib[0]);
    //добавляем пи к углу
    tetab[0] = -asin(sintb[0]) * 180 / M_PI + 180;

    for (int i = 0;i < n-1;i++) {
        //tb[i] = mb[i][0];
        ub[i+1] = mb[i][1];
        vb[i+1] = mb[i][2];
        rb[i+1] = mb[i][3];
        fib[i+1] = mb[i][4];
        pub[i+1] = mb[i][5];
        pvb[i+1] = mb[i][6];
        prb[i+1] = mb[i][7];
        pfib[i+1] = mb[i][8];
        hpb[i+1] = mb[i][9];
        sintb[i+1] = mb[i][10];
        costb[i+1] = mb[i][11];
        tb[i+1] = mb[i][12];
        tetab[i+1] = -asin(sintb[i+1])*180/M_PI+180;
        xb[i+1] = rb[i+1] * cos(fib[i+1]);
        yb[i+1] = rb[i+1] * sin(fib[i+1]);
    }

    for (int i = 0;i < n;i++)
        std::cout << tetab[i] << ";";

    // Prepare data.
    /*
    int n = 1000;
    std::vector<double> x(n), y(n), z(n), w(n, 2);
    for (int i = 0; i < n; ++i) {
        x.at(i) = i * i;
        y.at(i) = sin(2 * M_PI * i / 360.0);
        z.at(i) = log(i);
    }
    */
    /*
    std::cout << "tb = ";
    for (int i = 0;i < n;i++) {
        std::cout<<tb[i]<<";";
    }
    std::cout << std::endl;
    std::cout << "ub = ";
    for (int i = 0;i < n;i++) {
        std::cout << ub[i] << ";";
    }
    std::cout << std::endl;
    */

    // Set the size of output image = 1200x780 pixels
    //plt::figure_size(1200, 780);

    // Plot line from given x and y data. Color is selected automatically.
    plt::plot(tb, ub);

    // Plot a red dashed line from given x and y data.
    //plt::plot(x, w, "r--");

    // Plot a line whose name will show up as "log(x)" in the legend.
    //plt::named_plot("log(x)", x, z);

    // Set x-axis to interval [0,1000000]
    //plt::xlim(0, 1000 * 1000);

    // Add graph title
    plt::title("U(t)");

    // Enable legend.
    //plt::legend();

    // save figure
    const char* filename = "./u(t).png";
    plt::save(filename);

    plt::show();//clearing the plot

    plt::plot(tb, rb);
    plt::title("V(t)");
    filename = "./r(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, rb);
    plt::title("R(t)");
    filename = "./r(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, fib);
    plt::title("Fi(t)");
    filename = "./fi(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, pub);
    plt::title("PU(t)");
    filename = "./pu(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, pvb);
    plt::title("PV(t)");
    filename = "./pv(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, prb);
    plt::title("PR(t)");
    filename = "./pr(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, pfib);
    plt::title("PFi(t)");
    filename = "./pfi(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, hpb);
    plt::title("Hp(t)");
    filename = "./hp(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, sintb);
    plt::title("sin_teta(t)");
    filename = "./sin_teta(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, costb);
    plt::title("cos_teta(t)");
    filename = "./cos_teta(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, tetab);
    plt::title("teta(t)");
    filename = "./teta(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, xb);
    plt::title("x(t)");
    filename = "./x(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(tb, yb);
    plt::title("y(t)");
    filename = "./y(t).png";
    plt::save(filename);
    plt::show();

    plt::plot(xb, yb);
    plt::axis("equal");
    plt::xlim(-1*1e11, 1.5* 1e11);
    plt::ylim(-1 * 1e11, 1.5 * 1e11);
    plt::title("y(x)");
    filename = "./y(x).png";
    plt::save(filename);
    plt::show();
}
