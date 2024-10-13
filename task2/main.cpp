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
    return 1 - 4 * x - 8 * t - std::exp(x) + 0.027 * (24 * x * x + 2 + t * std::exp(x));
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
    const double h = 0.02;

    int n = (int)(1 / h) + 1;
    int j_max = (int) 1 / tau;

    matrix<double> u(j_max, n);

    for (int k = 0; k < n; k++){
        u.set_at(0, k, mu(k * h));
    }

    for (int j = 1; j < j_max; j++){
        u.set_at(j, 0, mu1(j * tau));
    }

    for (int j = 1; j < j_max; j++){
        u.set_at(j, n - 1, mu2(j * tau));
    }

    matrix<double> A(n - 1, n - 1);
    std::vector<double> b(n - 1);

    for (int i = 0; i < n - 1; i++){
        if (i < n - 2)
            A.set_at(i, i + 1, -a * tau);
        if (i > 0)
            A.set_at(i, i - 1, -a * tau);
        A.set_at(i, i, 2 * h * h + 2 * a * tau);
    }



    for (int j = 1; j < j_max; j++){
        for (int k = 0; k <= n - 3; k++){
            b[k] = u[j-1][k] * a * tau + u[j-1][k+1] * (2*h*h - 2*a*tau) + u[j-1][k+2]*a*tau + f((k + 1) * h, (j - 0.5)*tau)*tau*2*h*h;
        }
        b[0] += u[j][0]*a*tau;
        b[n-1] += u[j][n-1]*a*tau;


        std::vector<double> x = solve_system(A, b);

        for (int i = 0; i < n - 2 ; i++){
            u.set_at(j, i + 1, x[i]);
        }

    }
    
    for (int j = 0; j < j_max ; j++){
        for (int k = 0; k < n; k++){
            std::cout << k * h << '\t' << j * tau << '\t' << u[j][k] << '\n';
        }
    }
}