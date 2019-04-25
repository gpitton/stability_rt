#ifndef DARCY_HEADER
#define DARCY_HEADER
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
//#include <Spectra/GenEigsSolver.h>
#include "basis.h"
#include "utils.h"


// template parameters: basis function type, number of expansion modes
template<class T, int N>
class Darcy {
 public:
    Darcy(double,double,double,double,double,double=0);
    void precompute_matrices();
    void assemble_matrix();
    void solve_eigenproblem();
    inline void recompute_constants() { C = 3.*std::sqrt(2.)/4./Bo; };
    inline void set_biot(const double biot) { Bo = biot; };
    inline void set_peclet(const double peclet) { Pe = peclet; };
    inline int count_negative_eigenvalues();
    inline int count_positive_eigenvalues();
    Eigen::Array<std::complex<double>, N + (N - 2), 1>& get_eigenvalues();
    //std::complex<double>& get_maximum_eigenvalue();
 protected:
    void build_pressure_pressure_block();
    void build_pressure_phase_block();
    void build_phase_pressure_block();
    void build_phase_phase_block();
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
    // the following declaration is removed as the use of Li and Lj as member functions
    // would make shared-memory parallelism very hard
    // T         Li, Lj;
    // shape of phase field in the initial configuration
    Phi_0     phi_0;
    // chemical potential
    Eta       eta;
    // vector to store the eigenvalues
    Eigen::Array<std::complex<double>, N + Np, 1>   w;
    // A and B are respectively the lhs and the rhs for the generalized eigenvalue problem
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A, B;
    // parameter-dependent blocks of A, which are computed independently to speed-up families of parameter-dependent computations
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  A1_M,
                                   A1_D;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>   A2_sgn,
                                   A2_C;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>   A3;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>    A4_sgn,
                                   A4_C,
                                   A4_Pexx,
                                   A4_Peyy,
                                   M;
    //SchurOperator<N, N> schur;
    // solver for the generalized eigenvalue problem
    Eigen::GeneralizedEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> ges;
    // solver for the Schur complement eigenvalue problem
    //Spectra::GenEigsSolver<double, Spectra::SMALLEST_REAL, SchurOperator<N, N>> ses;
    //template<int M> friend class Trig3Darcy;
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
    //Li(1),
    //Lj(1),
    phi_0(std::sqrt(cahn)),
    eta(std::sqrt(cahn), mode_k, mode_l),
    A(N + (N - 2), N + (N - 2)),
    B(N + (N - 2), N + (N - 2)),
    A1_M((N - 2), (N - 2)),
    A1_D((N - 2), (N - 2)),
    A2_sgn((N - 2), N),
    A2_C((N - 2), N),
    A3(N, (N - 2)),
    A4_sgn(N, N),
    A4_C(N, N),
    A4_Pexx(N, N),
    A4_Peyy(N, N)
//ses(&schur, 3, 6)   // will compute just three eigenvalues, the third parameter must be greater than 2*nev
{
    sign = s;
    Bo = biot;
    Pe = peclet;
    Ch = cahn;
    k = mode_k;
    l = mode_l;
    eps = std::sqrt(Ch);

    // initialize B to zero
    B = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(N + Np, N + Np);
}


template<class T, int N>
void Darcy<T, N>::assemble_matrix() {
    // squared modulus of the component orthogonal to the gravity of the wavevector
    double kort = std::pow(k, 2) + std::pow(l, 2);

    // from Eigen's documentation:
    // since block is a templated member, the keyword template has to be used if the matrix type is also a template paramete
    A.template block<Np, Np>(0, 0) = kort*A1_M - A1_D;
    A.template block<Np, N>(0, Np) = sign/2.*A2_sgn + C/eps*A2_C;
    A.template block<N, Np>(Np, 0) = A3;
    A.template block<N, N>(Np, Np) = -sign/2.*A4_sgn - C/eps*A4_C + 1./Pe*(- kort*A4_Pexx + A4_Peyy);
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
    A1_M = B.template block<Np, Np>(Np, Np);
    #pragma omp parallel for
    for (int i = 0; i < Np; ++i) {
        //Li = T(i + 1);
        T Li {i + 1};
        for (int j = 0; j < Np; ++j) {
            //Lj = T(j + 1);
            T Lj {j + 1};

            auto f = [&Li, &Lj](const double x) { return Li(x)*Lj.deriv(2, x); };
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
    #pragma omp parallel for
    for (int i = 0; i < Np; ++i) {
        T Li {i + 1};
        for (int j = 0; j < N; ++j) {
            T Lj {j + 1};

            auto lj   = std::bind(Lj, std::placeholders::_1);
            auto ljd  = std::bind([&Lj](const double x) { return Lj.deriv(1, x); }, std::placeholders::_1);
            auto ljd2 = std::bind([&Lj](const double x) { return Lj.deriv(2, x); }, std::placeholders::_1);
            auto ljd3 = std::bind([&Lj](const double x) { return Lj.deriv(3, x); }, std::placeholders::_1);

            auto f = [this, &Li, &lj, &ljd, &ljd2, &ljd3](const double x) {
                return Li(x)*(phi_0.deriv(2, x)*eta(lj, ljd2, x) +
                              phi_0.deriv(1, x)*eta.deriv(lj, ljd, ljd3, x));
                               };
            A2_C(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);

            auto g = [&Li, &Lj](const double x) { return Li(x)*Lj.deriv(1, x); };
            A2_sgn(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, -1., 1., 5, 1.e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_phase_pressure_block() {
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        T Li {i + 1};
        for (int j = 0; j < Np; ++j) {
            T Lj {j + 1};

            auto f = [this, &Li, &Lj](const double x) { return Li(x)*phi_0.deriv(1, x)*Lj.deriv(1, x); };
            A3(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::build_phase_phase_block() {
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        T Li {i + 1};
        for (int j = 0; j < N; ++j) {
            T Lj {j + 1};

            auto lj = std::bind(Lj, std::placeholders::_1);
            auto ljd  = std::bind([&Lj](const double x) { return Lj.deriv(1, x); }, std::placeholders::_1);
            auto ljd2 = std::bind([&Lj](const double x) { return Lj.deriv(2, x); }, std::placeholders::_1);
            auto ljd4 = std::bind([&Lj](const double x) { return Lj.deriv(4, x); }, std::placeholders::_1);

            auto f = [this, &Li, &Lj](const double x) { return Li(x)*phi_0.deriv(1, x)*Lj(x); };
            A4_sgn(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);

            auto g = [this, &Li, lj, ljd2](const double x) {
                return Li(x)*std::pow(phi_0.deriv(1, x), 2)*eta(lj, ljd2, x);
                };
            A4_C(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, -1., 1., 5, 1.e-9);

            auto gx = [this, &Li, lj, ljd2](const double x) { return Li(x)*eta(lj, ljd2, x); };
            A4_Pexx(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(gx, -1., 1., 5, 1.e-9);

            auto h = [this, &Li, lj, ljd, ljd2, ljd4](const double x) { return Li(x)*eta.deriv_yy(lj, ljd, ljd2, ljd4, x); };
            A4_Peyy(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(h, -1., 1., 5, 1.e-9);
        }
    }
}


template<class T, int N>
void Darcy<T, N>::precompute_matrices() {
    // in this loop we fill B: B_ij = (Li, Lj)
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        T Li {i + 1};
        for (int j = 0; j <= i; ++j) {
            T Lj {j + 1};

            auto f = [&Li, &Lj](const double x) { return Li(x)*Lj(x); };
            B(Np + i, Np + j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);
        }
    }

    for (int i=0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            B(Np + i, Np + j) = B(Np + j, Np + i);

    build_pressure_pressure_block();
    build_pressure_phase_block();
    build_phase_pressure_block();
    build_phase_phase_block();
    std::cout << "matrix block computation complete." << std::endl;
}


template<class T, int N>
Eigen::Array<std::complex<double>, N + (N - 2), 1>& Darcy<T, N>::get_eigenvalues() {
    return w;
}


template<class T, int N>
    inline int Darcy<T, N>::count_negative_eigenvalues() {
    return std::count_if(&w(0), &w(w.rows() - 1), [](const std::complex<double> x) { return (x.real() < 0.); });
}


template<class T, int N>
inline int Darcy<T, N>::count_positive_eigenvalues() {
    return std::count_if(&w(0), &w(w.rows() - 1), [](const std::complex<double> x) { return (x.real() > 0.); });
}


// specialization for Trig3 basis functions
// this allows to skip a few computations, which we know are going to return zero
//see: https://stackoverflow.com/questions/36079774/invalid-use-of-incomplete-type-class-method-specialization
// https://stackoverflow.com/questions/30816813/invalid-use-of-incomplete-type-for-partial-template-specialization-c
template<int N>
class Trig3Darcy : public Darcy<Trig3, N> {
 public:
    Trig3Darcy(double,double,double,double,double,double=0);
    void precompute_matrices();
 private:
    void build_pressure_pressure_block();
    void build_pressure_phase_block();
    void build_phase_pressure_block();
    void build_phase_phase_block();
    static constexpr int Np = N - 2;    // number of pressure modes
};


template<int N>
Trig3Darcy<N>::Trig3Darcy(double s,
                          double biot,
                          double peclet,
                          double cahn,
                          double mode_k,
                          double mode_l)
:
Darcy<Trig3, N>::Darcy(s, biot, peclet, cahn, mode_k, mode_l)
{
    // initialize to zero a bunch of matrices
    Trig3Darcy::A1_M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(Np, Np);
    Trig3Darcy::A1_D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(Np, Np);
    Trig3Darcy::A2_sgn = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(Np, N);
    Trig3Darcy::A2_C = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(Np, N);
    Trig3Darcy::A3 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, Np);
    Trig3Darcy::A4_sgn = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, N);
    Trig3Darcy::A4_C = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, N);
    Trig3Darcy::A4_Pexx = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, N);
    Trig3Darcy::A4_Peyy = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, N);
}


template<int N>
void Trig3Darcy<N>::precompute_matrices() {
    // in this loop we fill B: B_ij = (Li, Lj)
    // for the Trig3 class, B is pentadiagonal, with the non-principal diagonals having the pattern: x-0-x-0-x-0-...
    // principal diagonal
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        Trig3 Li {i + 1};

        auto f = [&Li](const double x) { return std::pow(Li(x), 2); };
        Trig3Darcy::B(Np + i, Np + i) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);
    }
    // first outer diagonal
    #pragma omp parallel for
    for (int i = 0; i < N - 2; i += 2) {
        Trig3 Li {i + 1};
        Trig3 Lj {i + 3};

        auto f = [&Li, &Lj](const double x) { return Li(x)*Lj(x); };
        Trig3Darcy::B(Np + i, Np + i + 2) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);
    }
    // second outer diagonal
    #pragma omp parallel for
    for (int i = 0; i < N - 4; i += 2) {
        Trig3 Li {i + 1};
        Trig3 Lj {i + 5};

        auto f = [&Li, &Lj](const double x) { return Li(x)*Lj(x); };
        Trig3Darcy::B(Np + i, Np + i + 4) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);
    }

    for (int i=0; i < N; ++i)
        for (int j = 0; j < i; ++j)
            Trig3Darcy::B(Np + i, Np + j) = Trig3Darcy::B(Np + j, Np + i);

    build_pressure_pressure_block();
    build_pressure_phase_block();
    build_phase_pressure_block();
    build_phase_phase_block();
    std::cout << "matrix block computation complete." << std::endl;
}


template<int N>
void Trig3Darcy<N>::build_pressure_pressure_block() {
    Trig3Darcy::A1_M = Trig3Darcy::B.template block<Np, Np>(Np, Np);
    // A1_D is pentadiagonal as well (proof: compute <cos (n pi x), d^2 sin(n pi x) + ...>), and symmetric
    // main diagonal
    #pragma omp parallel for
    for (int i = 0; i < Np; ++i) {
        Trig3 Li {i + 1};

        auto f = [&Li](const double x) { return Li(x)*Li.deriv(2, x); };
        Trig3Darcy::A1_D(i, i) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1e-9);
        }
    // first outer diagonal
    #pragma omp parallel for
    for (int i = 0; i < Np - 2; i += 2) {
        Trig3 Li {i + 1};
        Trig3 Lj {i + 3};

        auto f = [&Li, &Lj](const double x) { return Li(x)*Lj.deriv(2, x); };
        Trig3Darcy::A1_D(i, i + 2) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1e-9);
    }
    // second outer diagonal
    #pragma omp parallel for
    for (int i = 0; i < Np - 4; i += 2) {
        Trig3 Li {i + 1};
        Trig3 Lj {i + 5};

        auto f = [&Li, &Lj](const double x) { return Li(x)*Lj.deriv(2, x); };
        Trig3Darcy::A1_D(i, i + 4) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1e-9);
    }

    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < i; ++j)
            Trig3Darcy::A1_D(i, j) = Trig3Darcy::A1_D(j, i);
}


template<int N>
void Trig3Darcy<N>::build_pressure_phase_block() {
    // A2_sgn is antisymmetric
    // A2_C has a checkerboard sparsity pattern (with 0 on the main diagonal)
    // TODO for N large, we are still computing a lot of zeros for A2_sgn
    #pragma omp parallel for
    for (int i = 0; i < Np; ++i) {
        int j_start = 0;
        if (i % 2 == 0)
            j_start = 1;
        Trig3 Li {i + 1};
        for (int j = j_start; j < N; j += 2) {
            Trig3 Lj {j + 1};

            auto lj   = std::bind(Lj, std::placeholders::_1);
            auto ljd  = std::bind([&Lj](const double x) { return Lj.deriv(1, x); }, std::placeholders::_1);
            auto ljd2 = std::bind([&Lj](const double x) { return Lj.deriv(2, x); }, std::placeholders::_1);
            auto ljd3 = std::bind([&Lj](const double x) { return Lj.deriv(3, x); }, std::placeholders::_1);

            auto f = [this, &Li, &lj, &ljd, &ljd2, &ljd3](const double x) {
                return Li(x)*(Trig3Darcy::phi_0.deriv(2, x)*Trig3Darcy::eta(lj, ljd2, x) +
                              Trig3Darcy::phi_0.deriv(1, x)*Trig3Darcy::eta.deriv(lj, ljd, ljd3, x));
                               };
            Trig3Darcy::A2_C(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);

            // exploit the antisymmetry of A2_sgn
            if (j < i)
                continue;

            auto g = [&Li, &Lj](const double x) { return Li(x)*Lj.deriv(1, x); };
            Trig3Darcy::A2_sgn(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, -1., 1., 5, 1.e-9);
        }
    }

    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < i; ++j)
            Trig3Darcy::A2_sgn(i, j) = -Trig3Darcy::A2_sgn(j, i);
}


template<int N>
void Trig3Darcy<N>::build_phase_pressure_block() {
    // A3 has a checkerboard pattern
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int j_start = 0;
        if (i % 2 == 0)
            j_start = 1;
        Trig3 Li {i + 1};
        for (int j = j_start; j < Np; j += 2) {
            Trig3 Lj {j + 1};

            auto f = [this, &Li, &Lj](const double x) { return Li(x)*Trig3Darcy::phi_0.deriv(1, x)*Lj.deriv(1, x); };
            Trig3Darcy::A3(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);
        }
    }
}


template<int N>
void Trig3Darcy<N>::build_phase_phase_block() {
    // matrices have a checkerboard pattern (with nonzeros on the main diagonal)
    // A4_sgn, A4_Pexx, and A4_Peyy are symmetric
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int j_start = 0;
        if (i % 2 == 1)
            j_start = 1;
        Trig3 Li {i + 1};
        for (int j = j_start; j < N; j += 2) {
            Trig3 Lj {j + 1};

            auto lj = std::bind(Lj, std::placeholders::_1);
            auto ljd  = std::bind([&Lj](const double x) { return Lj.deriv(1, x); }, std::placeholders::_1);
            auto ljd2 = std::bind([&Lj](const double x) { return Lj.deriv(2, x); }, std::placeholders::_1);
            auto ljd4 = std::bind([&Lj](const double x) { return Lj.deriv(4, x); }, std::placeholders::_1);

            auto g = [this, &Li, lj, ljd2](const double x) {
                return Li(x)*std::pow(Trig3Darcy::phi_0.deriv(1, x), 2)*Trig3Darcy::eta(lj, ljd2, x);
            };
            Trig3Darcy::A4_C(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(g, -1., 1., 5, 1.e-9);

            if (j > i)
                continue;

            auto f = [this, &Li, &Lj](const double x) { return Li(x)*Trig3Darcy::phi_0.deriv(1, x)*Lj(x); };
            Trig3Darcy::A4_sgn(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f, -1., 1., 5, 1.e-9);

            auto gx = [this, &Li, lj, ljd2](const double x) { return Li(x)*Trig3Darcy::eta(lj, ljd2, x); };
            Trig3Darcy::A4_Pexx(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(gx, -1., 1., 5, 1.e-9);

            auto h = [this, &Li, lj, ljd, ljd2, ljd4](const double x) { return Li(x)*Trig3Darcy::eta.deriv_yy(lj, ljd, ljd2, ljd4, x); };
            Trig3Darcy::A4_Peyy(i, j) = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(h, -1., 1., 5, 1.e-9);
        }
    }

    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            Trig3Darcy::A4_sgn(i, j) = Trig3Darcy::A4_sgn(j, i);

    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            Trig3Darcy::A4_Pexx(i, j) = Trig3Darcy::A4_Pexx(j, i);

    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            Trig3Darcy::A4_Peyy(i, j) = Trig3Darcy::A4_Peyy(j, i);
}


#endif
