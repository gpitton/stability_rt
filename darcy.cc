#include "darcy.h"

using namespace std::placeholders;
using namespace boost::math::quadrature;


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
    eta(std::sqrt(cahn), mode_k, mode_l)
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
    
}


template<class T, int N>
void Darcy<T, N>::solve_eigenproblem() {
    
}


template<class T, int N>
void Darcy<T, N>::build_pressure_pressure_block() {
    // pressure mass matrix: A1_M_{ij} = (Li, Lj)
    for (int i = 0; i < Np; ++i) {
        Li = Trig3(i + 1);
        for (int j = 0; j < Np; ++j) {
            Lj = Trig3(j + 1);

            auto f = [this](double x) { return Li(x)*Lj(x); };
            A1_M(i, j) = gauss_kronrod<double, 15>::integrate(
                           f,
                           -1.,    // integration interval limit, left
                           1.,     // integration interval limit, right
                           5,      // maximum depth of adaptive integration
                           1e-9    // relative tolerance
                           );
            f = [this](double x) { return Li(x)*Lj.deriv(2, x); };
            A1_D(i, j) = gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_pressure_phase_block() {
    for (int i = 0; i < Np; ++i) {
        Li = Trig3(i + 1);
        for (int j = 0; j < N; ++j) {
            Lj = Trig3(j + 1);

            auto lj   = std::bind(Lj, _1);
            auto ljd  = std::bind(Lj.deriv, 1, _1);
            auto ljd2 = std::bind(Lj.deriv, 2, _1);
            auto ljd3 = std::bind(Lj.deriv, 3, _1);

            auto f = [this, lj, ljd, ljd2, ljd3](double x) {
                return Li(x)*(phi_0.deriv(2, x)*eta(lj, ljd2, x) +
                              phi_0.deriv(1, x)*eta.deriv(lj, ljd, ljd3, x));
                               };
            A2_C(i, j) = gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);

            f = [this](double x) { return Li(x)*Lj.deriv(x); };
            A2_sgn(i, j) = gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_phase_pressure_block() {
    for (int i = 0; i < N; ++i) {
        Li = Trig3(i + 1);
        for (int j = 0; j < Np; ++j) {
            Lj = Trig3(j + 1);

            auto f = [this](double x) { return Li(x)*phi_0.deriv(1, x)*Lj.deriv(1, x); };
            A3(i, j) = gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_phase_phase_block() {
    for (int i = 0; i < N; ++i) {
        Li = Trig3(i + 1);
        for (int j = 0; j < N; ++j) {
            Lj = Trig3(j + 1);

            auto lj = std::bind(Lj, _1);

            auto f = [this](double x) { return Li(x)*phi_0.deriv(1, x)*Lj(x); };
            A4_sgn(i, j) = gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);

            f = [this, lj](double x) {
                return Li(x)*std::pow(phi_0.deriv(1, x), 2)*eta(lj, x);
                };
            A4_C(i, j) = gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);

            f = [this, lj](double x) { return Li(x)*eta(lj, x); };
            A4_Pexx(i, j) = gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);

            f = [this, lj](double x) { return Li(x)*eta.deriv_yy(lj, x); };
            A4_Peyy(i, j) = gauss_kronrod<double, 15>::integrate(f, -1, 1, 5, 1e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::precompute_matrices() {
    build_pressure_pressure_block();
    build_pressure_phase_block();
    build_phase_pressure_block();
    build_phase_phase_block();
    // fill B here B_ij = (Li, Lj)
}
