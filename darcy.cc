#include "darcy.h"


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
    
}
