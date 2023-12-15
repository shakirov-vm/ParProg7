#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <cmath>
#include <omp.h>

//#define CHUNK_SIZE (128)
const int N = 1023;
const int max_iter = 10000;
const int b_max = 10;
const double h = 1.0 / N;
const double mul_const = h * h / 12;

// g++ para.cpp -fopenmp

// https://algowiki-project.org/en/Methods_for_solving_tridiagonal_SLAEs#Cyclic_reduction_method
// a_i * x_i-1 + b_i * x_i + c_i * x_i+1 = f_i 
void diag_matr_solve(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, 
                     std::vector<double>& f, std::vector<double>& curr_y) {
    
#pragma omp parallel
{
    // Forward
    unsigned long long stride = 1;
    for (unsigned long long nn = N, low = 2; nn > 1; nn /= 2, low *= 2, stride *= 2) {

        #pragma omp for 
        for (int i = low - 1; i < N; i += stride * 2) {
            double alpha = -a[i] / b[i - stride];
            double gamma = -c[i] / b[i + stride];
            a[i] = alpha * a[i - stride];
            b[i] = alpha * c[i - stride] + b[i] + gamma * a[i + stride];
            c[i] = gamma * c[i + stride];
            f[i] = alpha * f[i - stride] + f[i] + gamma * f[i + stride];
        }
    }

#pragma omp barrier
    // Backward
    curr_y[N / 2] = f[N / 2] / b[N / 2];
    for (stride /= 2; stride >= 1; stride /= 2) {

        #pragma omp for 
        for (unsigned long long i = stride - 1; i < N; i += stride * 2) {
            curr_y[i] = f[i];
            if (i - stride > 0) curr_y[i] -= a[i] * curr_y[i - stride];
            if (i + stride < N) curr_y[i] -= c[i] * curr_y[i + stride];
            curr_y[i] /= b[i];
        }
    }
}

}

struct layout {
    
    std::vector<double>& x;
    std::vector<double>& a;
    std::vector<double>& b;
    std::vector<double>& c;
    std::vector<double>& f;

    layout(std::vector<double>& x_, std::vector<double>& a_, std::vector<double>& b_, 
            std::vector<double>& c_, std::vector<double>& f_) : 
            x(x_), a(a_), b(b_), c(c_), f(f_) {}
};

void calc_param(layout& lay, int curr_bound,
                std::vector<double>& prev_y, std::vector<double>& curr_y,
                std::vector<std::vector<double>>& solutions) {

    auto a = lay.a;
    auto b = lay.b;
    auto c = lay.c;
    auto x = lay.x;
    auto f = lay.f;

    for (int iters = 0; iters < max_iter; iters++) {

        a[0] = 0;
        b[0] = 1;
        c[0] = 0;

        f[0] = 1;

        #pragma omp parallel for 
        for (int k = 1; k < N - 1; k++) {

            double exp_k_minus_1, exp_k, exp_k_plus_1;
            
            exp_k_minus_1 = std::exp(-prev_y[k - 1]);
            exp_k = std::exp(-prev_y[k]);
            exp_k_plus_1  = std::exp(-prev_y[k + 1]);
        
            a[k] = 1 - mul_const * exp_k_plus_1;
            b[k] = -2 - 10.0 * mul_const * exp_k;
            c[k] = 1 - mul_const * exp_k_minus_1;

            f[k] = mul_const * (exp_k_plus_1 * (1 - prev_y[k + 1]) + 
                10 * exp_k * (1 - prev_y[k]) + exp_k_minus_1 * (1 - prev_y[k - 1]));
        }

        a[N - 1] = 0;
        b[N - 1] = 1;
        c[N - 1] = 0;

        f[N - 1] = curr_bound * 1.0 / b_max;
        
        diag_matr_solve(a, b, c, f, curr_y);

        prev_y = curr_y;
    }

    solutions[curr_bound] = curr_y;
}

void prepare_param(std::vector<double>& x, std::vector<double> prev_y) {
        
    for (int i = 0; i < N; i++) {
        x[i] = i * h;
        prev_y[i] = 0;
    }
}

int main(int argc, char **argv) {

    if (argc != 2) {
        printf("Need 2 arg\n");
        return 1;
    }

    int prog_num_threads = atoi(argv[1]);
    omp_set_num_threads(prog_num_threads);

    std::vector<double> x(N);
    std::vector<double> prev_y(N);
    std::vector<double> curr_y(N);
    std::vector<double> f(N);
    // matrix diagonals
    std::vector<double> a(N);
    std::vector<double> b(N);
    std::vector<double> c(N);

    std::vector<std::vector<double>> solutions(b_max);
    for (auto it : solutions) {
    	it = std::vector<double>(N);
    }

    struct layout lay(x, a, b, c, f);

    for (int curr_bound = 0; curr_bound < b_max; curr_bound++) {

        prepare_param(x, prev_y);

        auto&& start = std::chrono::high_resolution_clock::now();
        calc_param(lay, curr_bound, prev_y, curr_y, solutions);
        auto&& end = std::chrono::high_resolution_clock::now();

        auto&& passed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Time in ms: " << passed << ", b: " << curr_bound << std::endl;
    }

    std::ofstream out_file;
    out_file.open("result_para.txt");

    out_file << "x ";
    for (int curr_bound = 0; curr_bound < b_max; curr_bound++) {
            out_file << " " << "y_b_" << curr_bound;
    }    
    out_file << std::endl;

    for (int i = 0; i < N; i++) {
        out_file << x[i];
        for (int curr_bound = 0; curr_bound < b_max; curr_bound++) {
            out_file << " " << solutions[curr_bound][i];
        }
        out_file << std::endl;
    }
    out_file.close();
}