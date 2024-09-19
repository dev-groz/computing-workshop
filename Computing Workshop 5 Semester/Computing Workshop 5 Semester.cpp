#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

#include <iomanip>

const double tau = 0.1;
const double h = 0.2;

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

double explicit_scheme(double x_j, double t_n) {
    return (U(x_j, t_n + tau) - U(x_j, t_n)) / tau + a * (U(x_j, t_n) - U(x_j - h, t_n)) / h;
}

int main()
{
    std::cout << std::fixed << std::setprecision(10);

    for (double t_i = 0; t_i < 1; t_i += tau) {
        for (double x_i = 0; x_i < 1; x_i += h) {
            std::cout << x_i << "\t" << t_i << "\t" << abs(f(x_i, t_i) - explicit_scheme(x_i, t_i)) << "\n";

        }
    }

}
