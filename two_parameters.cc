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


// p1, p2 stand for "parameter 1" and "parameter 2"
// data contains the stability data parametrized by p1 and p2
template<typename T1, typename T2>
void output_stability_data(const T1& p1, const T1& p2, const T2& data, std::string fname) {
    std::ofstream output_file(fname);
    std::ostream_iterator<double> output_iterator(output_file, "\n");
    // first line contains the lengths of p1, p2, and data
    output_file << std::to_string(p1.size()) << " " << std::to_string(p2.size()) << " " << std::to_string(data.size()) << "\n";
    std::copy(p1.begin(), p1.end(), output_iterator);
    std::copy(p2.begin(), p2.end(), output_iterator);
    std::copy(data.begin(), data.end(), output_iterator);
    output_file.close();
    std::cout << "\n" << "written file " << fname << std::endl;
}


int main() {
    const int N = 64;
    std::string cahn = "0.1";
    double Ch = std::stod(cahn);
    int k = 1;

    // range of Biot number
    int n_bo = 100,
        n_pe = 100;

    auto ebo = linspace<double>(-1., 2., n_bo);
    auto Bo (ebo);
    std::transform(ebo.begin(), ebo.end(), Bo.begin(), [](double x) { return std::pow(10., x); });
    auto epe = linspace<double>(0., 4., n_pe);
    auto Pe (epe);
    std::transform(epe.begin(), epe.end(), Pe.begin(), [](double x) { return std::pow(10., x); });
    std::vector<int> res (n_bo*n_pe);


    Darcy<Trig3, N>     stability_problem(-1., Bo[0], Pe[0], Ch, k);
    stability_problem.precompute_matrices();

    for (int i = 0; i < Pe.size(); ++i) {
        stability_problem.set_peclet(Pe[i]);
        for (int j = 0; j < Bo.size(); ++j) {
            stability_problem.set_biot(Bo[j]);
            stability_problem.recompute_constants();
            stability_problem.assemble_matrix();
            stability_problem.solve_eigenproblem();
            res[i*n_bo + j] = N - stability_problem.count_negative_eigenvalues();
       }
    }

    output_stability_data<std::vector<double>, std::vector<int>> (Pe, Bo, res, "output_" + cahn + ".txt");

    return 0;
}

