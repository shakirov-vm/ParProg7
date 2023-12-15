#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <cmath>

const int N = 1023;
const int max_iter = 10000;
const int b_max = 10;
const double h = 1.0 / N;
const double mul_const = h * h / 12;
const double eps = 1e-10;

bool is_suff_accur(std::vector<double>& prev_y, std::vector<double>& curr_y) {
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += std::abs(prev_y[i] - curr_y[i]);
    }
    return res < eps;
}

void diag_matr_solve(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, 
                     std::vector<double>& f, std::vector<double>& x) {
    // Forward
    unsigned long long stride = 1;
    for (unsigned long long nn = N, low = 2; nn > 1; nn /= 2, low *= 2, stride *= 2) {
        for (int i = low - 1; i < N; i += stride * 2) {
            double alpha = -a[i] / b[i - stride];
            double gamma = -c[i] / b[i + stride];
            a[i] = alpha * a[i - stride];
            b[i] = alpha * c[i - stride] + b[i] + gamma * a[i + stride];
            c[i] = gamma * c[i + stride];
            f[i] = alpha * f[i - stride] + f[i] + gamma * f[i + stride];
        } 
    } 
    // Backward
    x[N / 2] = f[N / 2] / b[N / 2];
    for (stride /= 2; stride >= 1; stride /= 2) {
        for (unsigned long long i = stride - 1; i < N; i += stride * 2) {
            x[i] = f[i];
            if (i - stride > 0) x[i] -= a[i] * x[i - stride];
            if (i + stride < N) x[i] -= c[i] * x[i + stride];
            x[i] /= b[i];
        }
    }
}

int main() {

	std::vector<double> b_range(10);
	for (int i = 0; i < 10; i++) {
		b_range[i] = 1 / 10 * i;
	}

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

    for (int curr_bound = 0; curr_bound < b_max; curr_bound++) {

        // Initial approximation is: y = 1 + (b - 2) * x + x^2 
        for (int i = 0; i < N; i++) {
            x[i] = i * h;
            prev_y[i] = 1.0 + (curr_bound * 1.0 / b_max - 2.0) * x[i] + x[i] * x[i];
        }

        auto&& start = std::chrono::high_resolution_clock::now();
        for (int iters = 0; iters < max_iter; iters++) {

            a[0] = 0;
            b[0] = 1;
            c[0] = 0;

            f[0] = 1;

    		for (int k = 1; k < N - 1; k++) {

                double exp_k_minus_1 = std::exp(-prev_y[k - 1]);
                double exp_k = std::exp(-prev_y[k]);
                double exp_k_plus_1 = std::exp(-prev_y[k + 1]);
            
                a[k] = 1.0 - mul_const * exp_k_plus_1;
                b[k] = -2.0 - 10.0 * mul_const * exp_k;
                c[k] = 1.0 - mul_const * exp_k_minus_1;

                f[k] = mul_const * (exp_k_plus_1 * (1.0 - prev_y[k + 1]) + 
                	10.0 * exp_k * (1.0 - prev_y[k]) + exp_k_minus_1 * (1.0 - prev_y[k - 1]));
            }

            a[N - 1] = 0;
            b[N - 1] = 1;
            c[N - 1] = 0;

            f[N - 1] = curr_bound * 1.0 / b_max;
            
            diag_matr_solve(a, b, c, f, curr_y);

            if (is_suff_accur(prev_y, curr_y)) {
//                std::cout << "Converge earlier" << std::endl;
//                break;
            }
            prev_y = curr_y;
        }

        solutions[curr_bound] = curr_y;

        auto&& end = std::chrono::high_resolution_clock::now();
        auto&& passed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Time in ms: " << passed << ", b: " << curr_bound << std::endl;
    }
    // Just broke on 0.5
    std::ofstream out_file;
    out_file.open("result.txt");

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