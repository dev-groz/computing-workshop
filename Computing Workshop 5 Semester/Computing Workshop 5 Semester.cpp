#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <vector>

#include <iomanip>

double U(double x, double t) {
    return pow(x + 1.1 * t, 2) - sin(2 * M_PI * t) / 2 - 3.1 * t * x;
}

double f(double x_i, double t_i) {
    return 2.168 * t_i - 0.34 * x_i - M_PI * cos(2 * M_PI * t_i);
}

double phi(double x) {
    return x * x;
}

double g1(double t) {
    return 1.21 * t * t - sin(2 * M_PI * t) / 2;
}


// splot 'solution.dat', (x+1.1*y)**2 - sin(2*pi*y)/2 - 3.1*y*x

int main()
{
    const double tau = 0.01;
    const double h = 0.001;

    const double a = 0.28;

    std::ofstream solution_file("solution.dat");
    solution_file << std::fixed << std::setprecision(10);

    const int t_len = (int)(1 / tau);
    const int x_len = (int)(1 / h);


    double matrix[t_len][x_len] {0};

    double max_error = 0;

    for (int n = 0; n < t_len; ++n){
        for (int j = 0; j < x_len; ++j){
            matrix[n][j] = 0;
        }
    }

    for (int j = 0; j < x_len; ++j){
        matrix[0][j] = phi(h * j);
    }

    for (int n = 1; n < t_len; ++n){
        matrix[n][0] = g1(tau * n);
    }

    for (int n = 0; n < t_len - 1; ++n){
        for (int j = 1; j < x_len; ++j){
            double x_j = j * h;
            double t_n = n * tau;

            matrix[n + 1][j] = tau * (f(x_j, t_n) - a * (matrix[n][j] - matrix[n][j-1])/h) + matrix[n][j];
        }
    }



    for (int n = 0; n < t_len; ++n){
        for (int j = 0; j < x_len; ++j){
            double error = std::abs(U(j * h, n * tau) - matrix[n][j]);
            if (error > max_error){
                max_error = error;
            }
            solution_file << j * h << "\t" << n * tau << "\t" << matrix[n][j] << "\n";

        }
    }

    std::cout << "Max error: " << max_error << std::endl;

    solution_file.close();
}
