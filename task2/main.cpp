#include <iostream>
#include <cmath>
#include <vector>


class matrix{
    public:
    matrix(int height, int width): values(height, std::vector<double>(width, 0)){
        this->height = height;
        this->width = width;
    }
    
    std::vector<double> operator[](int index) const {
        return values[index];
    }

    void set_at(int y, int x, double value){
        values[y][x] = value;
    }

    int get_height(){
        return height;
    }

    int get_width(){
        return width;
    }

    private:
    int height;
    int width;
    std::vector<std::vector<double>> values;
};


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
    const double tau = 0.1;
    const double h = 0.05;

    int j_max = (int) 1 / tau;
    int k_max = (int) 1 / h;

    matrix grid(j_max, k_max);

    for (int k = 0; k < k_max; k++){
        grid.set_at(0, k, mu(k * h));
    }

    for (int j = 1; j < j_max; j++){
        grid.set_at(j, 0, mu1(j * tau));
    }

    for (int j = 1; j < j_max; j++){
        grid.set_at(j, k_max - 1, mu2(j * tau));
    }

    for (int j = 0; j < j_max ; j++){
        for (int k = 0; k < k_max; k++){
            std::cout << k * h << '\t' << j * tau << '\t' << grid[j][k] << '\n';
        }
    }
}