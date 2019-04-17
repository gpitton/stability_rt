#ifndef DARCY_HEADER
#define DARCY_HEADER
#include "basis.h"
#include "utils.h"


// template parameters: basis function type, number of expansion modes
template<class T, int N>
class Darcy {
 public:
    Darcy(double,double,double,double,double,double=0);
    void assemble_matrix();
    void solve_eigenproblem();
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
};

#endif
