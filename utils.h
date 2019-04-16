#ifndef UTILS_HEADER
#define UTILS_HEADER
#include <functional>


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

#endif
