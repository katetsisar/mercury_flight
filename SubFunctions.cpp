#include "SubFunctions.h"

// x - [pu,pv,pr,t1]

//fun2
inline void differentiate(int dimension,
    double x0,
    vector<double> y,
    vector<double>& dy,
    int& exit_signal,
    unsigned int step_number,
    Matrix& result) {
    //std::cout << "\nDIFF" << " step_num=" << step_number;
    //std::cout << "\nx0=" << x0;
    //вычисляем параметры на данный момент
    //std::cout << "y = [" << y[0] << " " << y[1]<<" ]";
    double sinTETA = y[4] / sqrt(y[4] * y[4] + y[5] * y[5]);
    double cosTETA = y[5] / sqrt(y[4] * y[4] + y[5] * y[5]);
    //double t1 = x0;
    double T = t1 * 86400;
    //std::cout << "\nsinTETA = " << sinTETA << "\ncosTETA = " << cosTETA << "\nt1 = " << t1<< "\nT = " << T;
    //производные
    //y=[u,v,r,fi,psi1,psi2,psi3,psi4]
    dy[0] = (y[1] * y[1] / y[2] - A0 * R0 * R0 / (y[2] * y[2]) + Totn * sinTETA / (1 - Qotn *x0* t1)) * T;
    dy[1] = (-y[0] * y[1] / y[2] + Totn * cosTETA / (1 - Qotn * x0*t1)) * T;
    dy[2] = y[0] * T;
    dy[3] = y[1] * T / y[2];
    dy[4] = (y[5] * y[1] / y[2] - y[6]) * T;
    dy[5] = (-2 * y[4] * y[1] / y[2] + y[5] * y[0] / y[2] - y[7] / y[2]) * T;
    dy[6] = (y[4] * y[1] * y[1] / (y[2]*y[2]) - 2 * y[4] * R0 * R0 * A0 / (y[2] * y[2] * y[2]) - y[5] * y[0] * y[1] / (y[2] * y[2]) + y[1] * y[7] / (y[2] * y[2])) * T;
    dy[7] = 0;
    //std::cout << "\ndiff = [" << dy[0] << " " << dy[1] << " " << dy[2] << " " << dy[3] << " " << dy[4] << " " << dy[5] << " " << dy[6] << " " << dy[7] << " ]";
    //вычисление Hpontr,tabs, result
    double Hpontr = (y[4]*dy[0]+y[5]*dy[1]+y[6]*dy[2]+y[7]*dy[3])/T - 1.0;
    //std::cout << "\nin diff " << Hpontr;
    tabs = x0 * t1;
    //std::cout << "\nHpontr= " << Hpontr << " tabs=" << tabs;
    //RES1
    //std::cout << "  dim = " << dimension;
    result[step_number + 1][dimension + 1] = Hpontr;//8+1
    result[step_number + 1][dimension + 2] = sinTETA;
    result[step_number + 1][dimension + 3] = cosTETA;
    result[step_number + 1][dimension + 4] = tabs; 
    //std::cout << "\nwriting result row = " << step_number + 1;
    //std::cout << "\n";
    //for (int i = 0; i < step_number + 1; i++)
    //{
        //for (int j = 0;j < dimension + 5;j++)
            //std::cout << result[i][j] << " ";
        //std::cout << "\n";
    //}
}


void getResidual(
    const vector<double>& x,
    vector<double>& residuals) {
    int exit_signal = 0;
    int step_counter = 0;
    int parametresNumber = 8;
    int steps = 1000;
    double firstTime = 0;
    double step = 1. / steps;
    //std::cout << "\nSTEP = "<<step;
    //std::cout << "\nSTEPS = " << steps;
    //std::cout << "called getRes\n";
    //double X0 = 0.0; == firstTime
    double PFi = 0.0;
    vector <double> x_res = { x[0] * (PU),x[1] * (PV),x[2] * (PR), x[3] * t1 };
    
    t1 = x[3] * 200;

    //std::cout << "\nx in residuals = [" << x_res[0] << " " << x_res[1] << " " << x_res[2]<<" ]\n";
    //t1 = x[3];
    double sinTETA = x_res[0] / sqrt(x_res[0] * x_res[0] + x_res[1] * x_res[1]);
    double cosTETA = x_res[1] / sqrt(x_res[0] * x_res[0] + x_res[1] * x_res[1]);
    double Hpontr = x_res[0] * (V0 * V0/R0-A0 + Totn * sinTETA) +
        x_res[1] * (-V0 * U0 / R0 + Totn * cosTETA ) +
        x_res[2] * U0  +
        PFi * (V0 / R0) - 1.0;//берем pfi вместо х[3], причем задаем знач = 0

    //double tabs = firstTime * x[3];
    double tabs = firstTime * t1;
    vector<double> y = { U0, V0, R0,FiZ,x_res[0],x_res[1],x_res[2],PFi };
    
    //запись RES1, создаем локальный
    Matrix result(steps + 1, vector<double>(parametresNumber + 5, 0));

    //result[0][dimension + 1] = Hpontr;
    //result[0][dimension + 2] = sinTETA;
    //result[0][dimension + 3] = cosTETA;
    //result[0][dimension + 4] = tabs;
    
    //вызов РК4
    RK4(differentiate, parametresNumber,steps, firstTime, step,y,exit_signal,step_counter, result);// в вызове передаем ссылку на функцию вычисления дифф

    //невязка F
    
    //changed y to result from rk4
    residuals[0] = result[steps][1];
    residuals[1] = (result[steps][2] - VM) / VM;
    residuals[2] = (result[steps][3] - RM) / RM;
    residuals[3] = result[steps][9];//Hpontr

}

double getEuclideanNorm(vector<double>& residuals){
    double norm = 0;
    for (int i = 0; i < residuals.size(); i++)
        norm += residuals[i] * residuals[i];
    norm = sqrt(norm);
    return norm;
}

double getLocalNorm(vector<double>& residuals, Matrix& A){
    double norm = 0;
    for (int i = 0; i < residuals.size(); i++)
    {
        double z = 0;
        for (int j = 0; j < residuals.size(); j++)
            z += A[i][j] * A[i][j];
        norm += residuals[i] * residuals[i] / z;
    }
    norm = sqrt(norm);
    return norm;
}



