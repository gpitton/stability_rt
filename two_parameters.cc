#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "darcy.h"


int main() {
    const int N = 64;
    std::string cahn = "0.1";
    double Ch = std::stod(cahn);
    int k = 1;

    // range of Biot number
    int n_bo = 100,
        n_pe = 100;

    // We distribute the parameters's range equally
    // on a log scale.
    auto ebo = linspace<double>(-1., 2., n_bo);
    auto Bo (ebo);
    std::transform(ebo.begin(), ebo.end(), Bo.begin(), [](double x) { return std::pow(10., x); });
    auto epe = linspace<double>(0., 4., n_pe);
    auto Pe (epe);
    std::transform(epe.begin(), epe.end(), Pe.begin(), [](double x) { return std::pow(10., x); });
    std::vector<int> res (n_bo*n_pe);


    //Darcy<Trig3, N>     stability_problem(-1., Bo[0], Pe[0], Ch, k);
    Trig3Darcy<N>     stability_problem(-1., Bo[0], Pe[0], Ch, k);
    stability_problem.precompute_matrices();

    for (int i = 0; i < Pe.size(); ++i) {
        stability_problem.set_peclet(Pe[i]);
        for (int j = 0; j < Bo.size(); ++j) {
            stability_problem.set_biot(Bo[j]);
            stability_problem.recompute_constants();
            stability_problem.assemble_matrix();
            stability_problem.solve_eigenproblem();
            res[i*n_bo + j] = stability_problem.count_positive_eigenvalues();
       }
    }

    output_stability_data<std::vector<double>, std::vector<int>> (Pe, Bo, res, "output_" + cahn + ".txt");

    return 0;
}

