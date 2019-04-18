#include <iostream>
#include "darcy.h"


int main() {
    const int N = 32;
    Darcy<Trig3, N>     stability_problem(1., 10., 1e3, 0.1, 1.);
    stability_problem.assemble_matrix();
    stability_problem.solve_eigenproblem();
    return 0;
}
