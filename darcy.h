#ifndef DARCY_HEADER
#define DARCY_HEADER
#include <cmath>
#include <complex>
#include <functional>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
//#include <Spectra/GenEigsSolver.h>
#include "basis.h"
#include "utils.h"


template<int R, int C>
class SchurOperator
{
 public:
    int rows() { return R; }
    int cols() { return C; }
    // y_out = M * x_in
    void perform_op(const double *x_in, double *y_out)
    {
        for(int i = 0; i < rows(); ++i)
            {
                y_out[i] = x_in[i] * (i + 1);
            }
    }
};


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
    // shape of phase field in the initial configuration
    Phi_0     phi_0;
    // chemical potential
    Eta       eta;
    // vector to store the eigenvalues
    Eigen::Array<std::complex<double>, N + Np, 1>   w;
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
    SchurOperator<N, N> schur;
    // solver for the generalized eigenvalue problem
    Eigen::GeneralizedEigenSolver<Eigen::Matrix<double, N + Np, N + Np>> ges;
    // solver for the Schur complement eigenvalue problem
    //Spectra::GenEigsSolver<double, Spectra::SMALLEST_REAL, SchurOperator<N, N>> ses;
};


// mode_l is initialized to zero by default: see constructor declaration
template<class T, int N>
Darcy<T, N>::Darcy(double s,
                   double biot,
                   double peclet,
                   double cahn,
                   double mode_k,
                   double mode_l)
    :
    phi_0(std::sqrt(cahn)),
    eta(std::sqrt(cahn), mode_k, mode_l),
    Li(1),
    Lj(1)
//ses(&schur, 3, 6)   // will compute just three eigenvalues, the third parameter must be greater than 2*nev
{
    sign = s;
    Bo = biot;
    Pe = peclet;
    Ch = cahn;
    k = mode_k;
    l = mode_l;
    eps = std::sqrt(Ch);
}


template<class T, int N>
void Darcy<T, N>::assemble_matrix() {
    // squared modulus of the component orthogonal to the gravity of the wavevector
    double kort = std::pow(k, 2) + std::pow(l, 2);

    A.block<Np, Np>(0, 0) = A1_M;
    //A.block<Np, Np>(0, 0) = kort*A1_M - A1_D;
    //A.block<Np, N>(0, Np) = sign/2.*A2_sgn + C/eps*A2_C;
    //A.block<N, Np>(Np, 0) = A3;
    //A.block<N, N>(Np, Np) = -sign/2.*A4_sgn - C/eps*A4_C + 1./Pe*(- kort*A4_Pexx + A4_Peyy);
}


template<class T, int N>
void Darcy<T, N>::solve_eigenproblem() {
    ges.compute(A, B);
    w = ges.eigenvalues().transpose();
}


template<class T, int N>
void Darcy<T, N>::build_pressure_pressure_block() {
    // pressure mass matrix: A1_M_{ij} = (Li, Lj)
    // this selects the sub-block of B of size <Np, Np>, starting at indices (Np-1, Np-1)
    A1_M = B.block<Np, Np>(Np - 1, Np - 1);
    for (int i = 0; i < Np; ++i) {
        Li = T(i + 1);
        for (int j = 0; j < Np; ++j) {
            Lj = T(j + 1);

            auto f = [this](double x) { return Li(x)*Lj.deriv(2, x); };
            A1_D(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(
                           f,
                           -1.,    // integration interval limit, left
                           1.,     // integration interval limit, right
                           5,      // maximum depth of adaptive integration
                           1e-9    // relative tolerance
                           );
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_pressure_phase_block() {
    for (int i = 0; i < Np; ++i) {
        Li = T(i + 1);
        for (int j = 0; j < N; ++j) {
            Lj = T(j + 1);

            auto lj   = std::bind(Lj, std::placeholders::_1);
            auto ljd  = std::bind(Lj.deriv, 1, std::placeholders::_1);
            auto ljd2 = std::bind(Lj.deriv, 2, std::placeholders::_1);
            auto ljd3 = std::bind(Lj.deriv, 3, std::placeholders::_1);

            auto f = [this, lj, ljd, ljd2, ljd3](double x) {
                return Li(x)*(phi_0.deriv(2, x)*eta(lj, ljd2, x) +
                              phi_0.deriv(1, x)*eta.deriv(lj, ljd, ljd3, x));
                               };
            A2_C(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);

            f = [this](double x) { return Li(x)*Lj.deriv(x); };
            A2_sgn(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_phase_pressure_block() {
    for (int i = 0; i < N; ++i) {
        Li = T(i + 1);
        for (int j = 0; j < Np; ++j) {
            Lj = T(j + 1);

            auto f = [this](double x) { return Li(x)*phi_0.deriv(1, x)*Lj.deriv(1, x); };
            A3(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_phase_phase_block() {
    for (int i = 0; i < N; ++i) {
        Li = T(i + 1);
        for (int j = 0; j < N; ++j) {
            Lj = T(j + 1);

            auto lj = std::bind(Lj, std::placeholders::_1);

            auto f = [this](double x) { return Li(x)*phi_0.deriv(1, x)*Lj(x); };
            A4_sgn(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);

            f = [this, lj](double x) {
                return Li(x)*std::pow(phi_0.deriv(1, x), 2)*eta(lj, x);
                };
            A4_C(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);

            f = [this, lj](double x) { return Li(x)*eta(lj, x); };
            A4_Pexx(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);

            f = [this, lj](double x) { return Li(x)*eta.deriv_yy(lj, x); };
            A4_Peyy(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::precompute_matrices() {
    double tmp;
    // in this loop we fill B: B_ij = (Li, Lj)
    for (int i=0; i < N; ++i) {
        Li = T(i + 1);
        for (int j = 0; j < N; ++j) {
            if (j <= i) {
                Lj = T(j + 1);

                auto f = [this](double x) { return Li(x)*Lj(x); };
                tmp = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);
            } else {
                tmp = B(j, i);
            }
            B(i, j) = tmp;
        }
    }
    build_pressure_pressure_block();
    build_pressure_phase_block();
    build_phase_pressure_block();
    build_phase_phase_block();
}


#endif
