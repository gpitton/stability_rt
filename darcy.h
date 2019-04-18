#ifndef DARCY_HEADER
#define DARCY_HEADER
#include <cmath>
#include <Eigen/Dense>
#include "basis.h"
#include "utils.h"


// template parameters: basis function type, number of expansion modes
template<class T, int N>
class Darcy {
 public:
    Darcy(double,double,double,double,double,double=0);
    void assemble_matrix();
    void solve_eigenproblem();
    inline void recompute_constants() { C = 3.*std::sqrt(2.)/4./Bo; };
    inline void set_biot(const double biot) { Bo = biot; };
    inline void set_peclet(const double peclet) { Pe = peclet; };
 private:
    void build_pressure_pressure_block();
    void build_pressure_phase_block();
    void build_phase_pressure_block();
    void build_phase_phase_block();
    void precompute_matrices();
    double    sign,   // gravity direction: +1 or -1
              Bo,     // Biot number
              Pe,     // Peclet number
              Ch,     // Cahn number
              k,      // x-component of perturbation wavevector
              l,      // z-component of perturbation wavevector (= 0 for 2d problems)
              eps;    // = sqrt(Ch), interface width parameter
    double    C;      // constant appearing in the term eta*grad(phi) in the momentum equation
    static constexpr int Np = N - 2;    // number of pressure modes
    // choice of basis functions for the perturbation expansions
    // Li is understood to be the test function, Lj to be the trial function
    T         Li, Lj;
    // A and B are respectively the lhs and the rhs for the generalized eigenvalue problem
    Eigen::Matrix<double, N + Np, N + Np> A, B;
    // parameter-dependent blocks of A, which are computed independently to speed-up families of parameter-dependent computations
    Eigen::Matrix<double, Np, Np>  A1_M,
                                   A1_D;
    Eigen::Matrix<double, Np, N>   A2_sgn,
                                   A2_C;
    Eigen::Matrix<double, N, Np>   A3;
    Eigen::Matrix<double, N, N>    A4_sgn,
                                   A4_C,
                                   A4_Pexx,
                                   A4_Peyy,
                                   M;
};

#endif
