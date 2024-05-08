#include "Modyfied_Newton.h"

void newtonModyfied(
    void (*getResidual)(
        const vector<double> & x,
        vector<double> & residuals),
    unsigned int dimension,
    vector<double> & x,
    vector<double>& residuals, int& iteration) 
{
    //точность
    double E = 1e-06;
    int max_iterations = 100;
    iteration = 0;

    for (; iteration < max_iterations;) 
	{
		//посчитали невязку (по исходным данным)

        getResidual(x, residuals);

        // Проверка сходимости (в первый раз скорее всего 
		// не пройдёт, вряд ли попали в нужное решение)
        int R1 = 0;
        for (int i = 0; i < dimension; i++)
            if (fabs(residuals[i]) > E)
                R1 = 1;

        if (R1 == 0)
            return;
        // конец проверки сходимости

        vector<double> B = residuals;
        for (int i = 0; i < B.size(); i++) 
            B[i] = -B[i];

        // Вычисление матрицы Якоби
        Matrix A(dimension, vector<double>(dimension, 0));
        vector<double> x_temp = x;
        for (int i = 0; i < dimension; i++) 
		{
            double jStep = E * fabs(x_temp[i]); // шаг по i-ой переменной
            if (fabs(jStep) < 1e-10)
                jStep = E;
            x[i] = x_temp[i] + jStep;
            getResidual(x, residuals);
            // цикл по вычислению производных по невязкам
            for (int j = 0; j < dimension; j++)
                A[j][i] = (residuals[j] + B[j]) / jStep;
            x[i] = x_temp[i];
        }

        iteration++;

        // решение СЛАУ
        Matrix LU = A;
        vector<int> T;
        string error = "";

        LUP(A, LU, T, error);
        if (!error.empty()) {
            std::cout << "systems matrix is special" << std::endl;
            //break;
        }


        vector<double> dx(B.size(), 0);
        solveSystem(LU, B, T, dx);


        R1 = 0;

        for (int i = 0; i < dimension; i++) 
		{
			int z;
			if (fabs(x[i]) < E)
				z = 1;
			else
				z = x[i];
			// проверка изменения аргументов на Ньютоновском шаге на величину, 
			// соизмеримую с заданной точностью
			if (fabs(dx[i] / z) > E)
				R1 = 1;
        }
		if (R1 != 1)
			return;

		vector<double> euclideanNorm(100);
		vector<double> localNorm(100);
		for (int i = 0; i < 100; i++)
		{
			for (int j = 0; j < dimension; j++)
				x[j] += dx[j] / 100;
			getResidual(x, residuals);
			euclideanNorm[i] = getEuclideanNorm(residuals);
			localNorm[i] = getLocalNorm(residuals, A);
			if (i > 0)
				if (localNorm[i] > localNorm[i - 1]) //min locnorm
				{
					for (int j = 0; j < dimension; j++)
						x[j] -= dx[j] / 100;
					break;
				}
		}
        //std::cout << "\niter=" << iteration;
    }

}