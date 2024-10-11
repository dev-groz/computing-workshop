#include <iostream>
#include <cmath>
#include <vector>


template <typename T>
class matrix{
    public:
    matrix(int height, int width): values(height, std::vector<T>(width, 0)){
        this->height = height;
        this->width = width;
    }
    
    std::vector<T> operator[](int index) const {
        return values[index];
    }

    void set_at(int y, int x, T value){
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
    std::vector<std::vector<T>> values;
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

std::vector<double> solve_system(const matrix<double> a, const std::vector<double>& b){
    int n = b.size();

    std::vector<double> x(n, 0);
    
    std::vector<double> v(n, 0);
    std::vector<double> u(n, 0);

    v[0] = a[0][1] / (-a[0][0]);
    u[0] = (-b[0]) / (-a[0][0]);

    for (int i = 1; i < n - 1; i++){
        v[i] = a[i][i+1] / (-a[i][i] - a[i][i-1]*v[i-1] );
        u[i] = (a[i][i-1] * u[i-1] - b[i])/ (-a[i][i] - a[i][i-1]*v[i-1] );
    }

    v[n - 1] = 0;
    u[n - 1] = (a[n-1][n-2] * u[n-2] - b[n-1]) / (-a[n-1][n-1] - a[n-1][n-2]*v[n-2]);


    x[n-1] = u[n-1];

    for (int i = n-1; i > 0; i--){
        x[i-1] = v[i - 1] * x[i] + u[i-1];
    }

    return x;
}

int main(){
    const double a = 0.027;
    const double tau = 0.1;
    const double h = 0.05;

    int j_max = (int) 1 / tau;
    int k_max = (int) 1 / h;


    matrix<double> grid(j_max, k_max);

    for (int k = 0; k < k_max; k++){
        grid.set_at(0, k, mu(k * h));
    }

    for (int j = 1; j < j_max; j++){
        grid.set_at(j, 0, mu1(j * tau));
    }

    for (int j = 1; j < j_max; j++){
        grid.set_at(j, k_max - 1, mu2(j * tau));
    }


    int n = k_max - 2;

    matrix<double> A(n, n);
    std::vector<double> b(n);

    for (int i = 0; i < n; i++){
        if (i < k_max - 3)
            A.set_at(i, i + 1, -a * tau);
        if (i > 0)
            A.set_at(i, i - 1, -a * tau);
        A.set_at(i, i, 2 * h * h - 2 * a * tau);
    }

    int j = 1;

    for (int k = 2; k < n; k++){
        b[k - 1] = grid[j][k - 1] * a * tau + grid[j][k] * (2*h*h - a * tau) + grid[j][k + 1] + f(h * k, (j + 0.5) * tau) * tau * 2*h*h;
    }

    b[0] = grid[j][0] * a * tau + grid[j][1] * (2*h*h - a * tau) + grid[j][2] + f(h * 1, (j + 0.5) * tau) * tau * 2*h*h + grid[j + 1][0] * a * tau;
    b[n - 1] = grid[j][n - 2] * a * tau + grid[j][n - 1] * (2*h*h - a * tau) + grid[j][n] + f(h * (n - 1), (j + 0.5) * tau) * tau * 2*h*h + grid[j + 1][n] * a * tau;

    std::vector<double> x = solve_system(A, b);


    for (int i = 0; i < n ; i++){
        grid.set_at(j, i + 1, x[i]);
    }


    // for (int i = 0; i < k_max - 2 ; i++){
    //     for (int j = 0; j < k_max - 2; j++){
    //         std::cout << A[i][j] << ' ';
    //     }
    //     std::cout << "\n";
    // }
    

    
    for (int j = 0; j < j_max ; j++){
        for (int k = 0; k < k_max; k++){
            std::cout << k * h << '\t' << j * tau << '\t' << grid[j][k] << '\n';
        }
    }

}