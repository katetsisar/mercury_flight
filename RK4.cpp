#include "RK4.h"

void RK4(
    void (*differentiate)(
        int dimension,
        double x0,
        vector<double> y,
        vector<double>& dy,
        int& exit_signal,
        unsigned int step_number,
		Matrix& result),
    unsigned int dimension,
    int steps,
    double x0,
    double dx,
    vector<double> y,
    int& exit_signal,
    int& step_counter,
    Matrix& result)
{
    step_counter = 0;
    double x = x0;
    exit_signal = 0;

    vector<double> ak1(dimension, 0);
    vector<double> ak2(dimension, 0);
    vector<double> ak3(dimension, 0);
    vector<double> ak4(dimension, 0);
    vector<double> y1(dimension, 0);

    result = Matrix(steps + 1, vector<double>(dimension + 5, 0));
    result[0][0] = x;
    for (int i = 0; i < dimension; i++)
        result[0][i + 1] = y[i];

    vector<double> dy(dimension, 0);

    for (int i = 0; i < steps; i++)
	{
        differentiate(dimension, x, y, dy, exit_signal, step_counter, result);
        if (exit_signal == 1)
            break;
        for (int j = 0; j < dimension; j++) 
		{
            ak1[j] = dx * dy[j];
            y1[j] = y[j] + ak1[j] / 2.0;
        }
        double x1 = x + dx / 2.0;
        differentiate(dimension, x1, y1, dy, exit_signal, step_counter, result);
        for (int j = 0; j < dimension; j++) 
		{
            ak2[j] = dx * dy[j];
            y1[j] = y[j] + ak2[j] / 2.0;
        }

        differentiate(dimension, x1, y1, dy, exit_signal, step_counter, result);
        for (int j = 0; j < dimension; j++) 
		{
            ak3[j] = dx * dy[j];
            y1[j] = y[j] + ak3[j];
        }

        differentiate(dimension, x + dx, y1, dy, exit_signal, step_counter, result);

        result[i + 1][0] = x + dx;

		for (int j = 0; j < dimension; j++) 
		{
            ak4[j] = dx * dy[j];
            y[j] += (ak1[j] + 2.0 * ak2[j] + 2.0 * ak3[j] + ak4[j]) / 6.0;
            result[i + 1][j + 1] = y[j];
        }

        x += dx;
        step_counter++;
    }

    std::ofstream myfile2;
    myfile2.open("rk4_new.csv");
    for (int i = 0;i < steps + 1;i++) {
        for (int j = 0;j < 8 + 5;j++)
           myfile2<< result[i][j] << ";";
        myfile2 << "\n";
    }
    
}