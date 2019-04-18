#include "darcy.h"

using namespace boost::math::quadrature;


// mode_l is initialized to zero by default: see constructor declaration
template<class T, int N>
Darcy<T, N>::Darcy(double s,
                   double biot,
                   double peclet,
                   double cahn,
                   double mode_k,
                   double mode_l) {
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
    for (int i = 0; i < Np; ++i)
        A1_M(i, i) = 1.;
    for (int i = 0; i < Np; ++i) {
        Li = Trig3(i + 1);
        for (int j = 0; j < Np; ++j) {
            Lj = Trig3(j + 1);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_pressure_phase_block() {
    
}


template<class T, int N>
void Darcy<T, N>::build_phase_pressure_block() {
    
}


template<class T, int N>
void Darcy<T, N>::build_phase_phase_block() {
    
}


template<class T, int N>
void Darcy<T, N>::precompute_matrices() {
    build_pressure_pressure_block();
    build_pressure_phase_block();
    build_phase_pressure_block();
    build_phase_phase_block();
    // fill B here B_ij = (Li, Lj)
}
