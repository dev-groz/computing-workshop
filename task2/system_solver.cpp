#include <iostream>
#include <vector>
#include <cmath>

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



void print_mat(matrix a, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            std::cout << a[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

void print_vec(std::vector<double>& x){
    for (int i = 0; i < x.size(); i++){
        std::cout << x[i] << '\n';
    }
}

std::vector<double> solve_system(const matrix a, const std::vector<double>& b){
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
    const int n = 5;
    matrix A(n, n);
    std::vector<double> b(n, 0);

    for (int i = 0; i < n; i++){
        if (i < n - 1)
            A.set_at(i, i + 1, 1);
        if (i > 0)
            A.set_at(i, i - 1, 1);
        A.set_at(i, i, 3);
    }

    for (int i = 0; i < n; i++){
        b[i] = 1;
    }


    print_mat(A, n);
    print_vec(b);

    std::vector<double> x = solve_system(A, b);

    print_vec(x);
}