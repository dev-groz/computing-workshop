#include <iostream>
#include <cmath>

double mu(double x){
    return -2 * std::pow(x, 4) - x * x;
}

double mu1(double t){
    return - 4 * t * t;
}

double mu2(double t){
    return t * (1 - std::exp(1)) - std::pow(1 + 2 * t, 2) - 2;
}

double f(double x, double t){
    return 1.054 + x * (0.648 * x - 4) + std::exp(x) * (0.027 * t - 1);
}


int main(){
    const double a = 0.027;
    const double tau = 0.01;
    const double h = 0.02;

    std::cout << "a = " << a << '\n';
    std::cout << "tau = " << tau << '\n';
    std::cout << "h = " << h << '\n';
}