#ifndef UTILS_HEADER
#define UTILS_HEADER
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>


// equally spaced points betwen a and b, inspired from
// numpy's omonimous function
template<typename T>
std::vector<T> linspace(T a, T b, int n) {
    std::vector<T> dst(n);
    T h = (b - a)/static_cast<T> (n - 1);
    for (int i = 0; i < n; ++i)
        dst[i] = a + h*i;
    return dst;
}


// p1, p2 stand for "parameter 1" and "parameter 2"
// data contains the stability data parametrized by p1 and p2
template<typename T1, typename T2>
void output_stability_data(const T1& p1, const T1& p2, const T2& data, std::string fname) {
    std::ofstream output_file(fname);
    std::ostream_iterator<double> output_iterator(output_file, "\n");
    // first line contains the lengths of p1, p2, and data
    output_file << std::to_string(p1.size()) << " " << std::to_string(p2.size()) << " " << std::to_string(data.size()) << "\n";
    std::copy(p1.begin(), p1.end(), output_iterator);
    std::copy(p2.begin(), p2.end(), output_iterator);
    std::copy(data.begin(), data.end(), output_iterator);
    output_file.close();
    std::cout << "\n" << "written file " << fname << std::endl;
}


class Phi_0 {
public:
    Phi_0(double);
    double operator()(const double);
    double deriv(int, const double);
private:
    double eps;
};


class Eta {
 public:
    Eta(double, int, int);
    double operator()(std::function<double (const double)>,
                      std::function<double (const double)>,
                      const double);
    double deriv(std::function<double (const double)>,
                 std::function<double (const double)>,
                 std::function<double (const double)>,
                 const double);
    double deriv_xx(std::function<double (const double)>,
                    std::function<double (const double)>,
                    const double);
    double deriv_yy(std::function<double (const double)>,
                    std::function<double (const double)>,
                    std::function<double (const double)>,
                    std::function<double (const double)>,
                    const double);
 private:
    double eps;
    int    k, l;
    Phi_0  phi_0;
};


// Phi_0 describes the profile for the phase field in the initial configuration
Phi_0::Phi_0(double epsilon){
    eps = epsilon;
}


inline double Phi_0::operator()(const double x) {
    return std::tanh(x/std::sqrt(2.)/eps);
}

// derivatives of the phase field in the initial configuration
double Phi_0::deriv(int m, const double x) {
    if (m == 2)
        return -std::tanh(x/std::sqrt(2.)/eps)/std::pow(eps*std::cosh(x/std::sqrt(2.)/eps), 2);
    return 1./(std::sqrt(2.)*eps*std::pow(std::cosh(x/std::sqrt(2.)/eps), 2));
}


// Eta describes the chemical potential
Eta::Eta(double epsilon, int kx, int kz=0)
    :
    phi_0(epsilon)
{
    eps = epsilon;
    // x-component of the wave vector of the perturbation
    k = kx;
    // z-component of the wave vector of the perturbation (for a 2-d problem, this is zero)
    l = kz;
}


// evaluate the Frechet differential of Eta at the point x and with respect to the variation f
//     fd2 is the second derivative of the variation f
inline double Eta::operator()(std::function<double (const double)> f,
                              std::function<double (const double)> fd2,
                              const double x) {
    return 3*std::pow(phi_0(x), 2)*f(x) - (1. - std::pow(eps*k, 2) - std::pow(eps*l, 2))*f(x) - std::pow(eps, 2)*fd2(x);
}


// evaluate the derivatives of the Frechet differential of Eta at the point x with respect to the variation f.
//     fd is the derivative of the variation f
//     fd3 is the third derivative of the variation f
double Eta::deriv(std::function<double (const double)> f,
                  std::function<double (const double)> fd,
                  std::function<double (const double)> fd3,
                  const double x) {
    double p  = phi_0(x);
    double py = phi_0.deriv(1, x);
    return 6*p*py*f(x) + 3*std::pow(p, 2)*fd(x) - (1. - std::pow(eps*k, 2) - std::pow(eps*l, 2))*fd(x) - std::pow(eps, 2)*fd3(x);
}


// in 3-d, deriv_xx is the sum of the second derivatives in the two directions orthogonal to gravity.
//     f is the direction along which to evaluate the Frechet differential of Eta
//     fd2 is the second derivative of the direction along which to evaluate the Frechet differential
double Eta::deriv_xx(std::function<double (const double)> f,
                     std::function<double (const double)> fd2,
                     const double x) {
    return -(std::pow(k, 2) + std::pow(l, 2))*(3*std::pow(phi_0(x), 2)*f(x) - (1. - std::pow(eps*k, 2) - std::pow(eps*l, 2))*f(x) - std::pow(eps, 2)*fd2(x));
}


// second derivative of Eta in the gravity direction.
//     f is the direction for the Frechet differential
//     fd, fd2, fd4 are the derivatives of order 1, 2, and 4 of the the direction along which the Frechet differential of Eta is evaluated
double Eta::deriv_yy(std::function<double (const double)> f,
                     std::function<double (const double)> fd,
                     std::function<double (const double)> fd2,
                     std::function<double (const double)> fd4,
                     const double x) {
    double p   = phi_0(x);
    double py  = phi_0.deriv(1, x);
    double pyy = phi_0.deriv(2, x);
    return 6*(std::pow(py, 2) + p*pyy)*f(x) + 6*p*py*fd(x) + (3*std::pow(p, 2) - 1.)*fd2(x) - std::pow(eps, 2)*(-(std::pow(k, 2) + std::pow(l, 2))*fd2(x) + fd4(x));
}


#endif
