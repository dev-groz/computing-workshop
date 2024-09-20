#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>

#include <iomanip>

const double tau = 0.05;
const double h = 0.05;

const double a = 0.28;

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
    return 1.21 * t - sin(2 * M_PI * t) / 2;
}



int main()
{

    std::ofstream solution_file("solution.dat");
    solution_file << std::fixed << std::setprecision(10);

    int x_len = (int)(1/tau);
    int t_len = (int)(1/h);

    double matrix[t_len][x_len] {0};

    double max_error = 0;

    for (int n = 0; n < t_len; ++n){
        for (int j = 0; j < x_len; ++j){
            matrix[n][j] = 0;
        }
    }

    for (int n = 0; n < t_len; ++n){
        matrix[n][0] = phi(h * n);
    }

    for (int j = 0; j < x_len; ++j){
        matrix[0][j] = g1(tau * j);
    }

    for (int n = 0; n < t_len - 1; ++n){
        for (int j = 1; j < x_len; ++j){
            double x_j = j * h;
            double t_n = n * tau;
            matrix[n + 1][j] = tau * (f(x_j, t_n) - a * (matrix[n][j] - matrix[n][j-1])/h) + matrix[n][j];
        }
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
            max_error = std::max(max_error, std::abs(U(j * h, n * tau) - matrix[n][j]));
            solution_file << j * h << "\t" << n * tau << "\t" << matrix[n][j] << "\n";

        }
    }

    std::cout << "Max error: " << max_error << std::endl;

    solution_file.close();
}
