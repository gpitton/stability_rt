#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "darcy.h"


template<typename T>
std::vector<T> linspace(T a, T b, int n) {
    std::vector<T> dst(n);
    T h = (b - a)/static_cast<T> (n - 1);
    for (int i = 0; i < n; ++i)
        dst[i] = a + h*i;
    return dst;
}


int main() {
    const int N = 64;
    double Ch = 0.01;
    int k = 1;

    // range of Biot number
    int n_bo = 100,
        n_pe = 100;

    auto ebo = linspace<double>(1.e-2, 2., n_bo);
    auto Bo (ebo);
    std::transform(ebo.begin(), ebo.end(), Bo.begin(), [](double x) { return std::pow(10., x); });
    auto epe = linspace<double>(0., 4., n_pe);
    auto Pe (epe);
    std::transform(epe.begin(), epe.end(), Pe.begin(), [](double x) { return std::pow(10., x); });
    std::vector<int> res (n_bo*n_pe);


    Darcy<Trig3, N>     stability_problem(1., Bo[0], Pe[0], Ch, k);
    stability_problem.precompute_matrices();

    for (int i = 0; i < Pe.size(); ++i) {
        stability_problem.set_peclet(Pe[i]);
        for (int j = 0; j < Bo.size(); ++j) {
            std::cout << i << " " << j << std::endl;
            stability_problem.set_biot(Bo[j]);
            stability_problem.recompute_constants();
            stability_problem.assemble_matrix();
            stability_problem.solve_eigenproblem();
            res[i*n_bo + j] = stability_problem.count_negative_eigenvalues();
            std::cout << res[i*n_bo + j] << " ";
       }
    }

    std::ofstream output_file("res_" + std::to_string(Ch) + ".txt");
    std::ostream_iterator<int> output_iterator(output_file, "\n");
    std::copy(res.begin(), res.end(), output_iterator);

    return 0;
}

