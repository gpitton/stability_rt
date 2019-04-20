#include <iostream>
#include "darcy.h"


int main() {
    const int N = 64;
    Darcy<Trig3, N>     stability_problem(1., 10., 1e3, 0.1, 1.);
    stability_problem.precompute_matrices();
    stability_problem.assemble_matrix();
    stability_problem.solve_eigenproblem();
    std::cout << "eigenvalues: " << stability_problem.get_maximum_eigenvalue() << std::endl;

    Darcy<Trig3, 2*N>     st(1., 10., 1e3, 0.1, 1.);
    st.precompute_matrices();
    st.assemble_matrix();
    st.solve_eigenproblem();
    std::cout << "eigenvalues: " << st.get_maximum_eigenvalue() << std::endl;

    return 0;
}
