#include<cmath>
#include "basis.h"


// class for defining a set of orthogonal basis functions on [-1, 1]
// such that their first and third derivatives are zero at x = \pm 1
Trig3::Trig3(int num) {
    if (num % 2 == 0) {
        n = num/2;
        t = Even;
    } else {
        n = (num + 1)/2;
        t = Odd;
        an = 4.*num/(2.*num + 3.);
        bn = num*(2.*num + 1.)/((num + 2.)*(2.*num + 3.));
    }
}


double Trig3::operator()(const double x) {
    if (t == Even) {
        return std::cos(n*pi*x);
    } else {
        return std::sin(n*pi*x) + an*std::sin((n + 1)*pi*x) + bn*std::sin((n + 2)*pi*x);
    }
}


double Trig3::deriv(int m=1, const double x=0.) {
    if (t == Even) {
        switch (m) {
        case 1:
            return -n*pi*std::sin(n*pi*x);
        case 2:
            return -std::pow(n*pi, 2)*std::cos(n*pi*x);
        case 3:
            return std::pow(n*pi, 3)*std::sin(n*pi*x);
        case 4:
            return std::pow(n*pi, 4)*std::cos(n*pi*x);
        }
    } else {
        switch (m) {
        case 1:
            return n*pi*(std::cos(n*pi*x) + (n + 1.)/n*an*std::cos((n + 1)*pi*x) + (n + 2.)/n*bn*std::cos((n + 2)*pi*x));
        case 2:
            return -std::pow(n*pi, 2)*(std::sin(n*pi*x) + std::pow((n + 1.)/n, 2)*an*std::sin((n + 1)*pi*x) + std::pow((n + 2.)/n, 2)*bn*std::sin((n + 2)*pi*x));
        case 3:
            return -std::pow(n*pi, 3)*(std::cos(n*pi*x) + std::pow((n + 1.)/n, 3)*an*std::cos((n + 1)*pi*x) + std::pow((n + 2.)/n, 3)*bn*std::cos((n + 2)*pi*x));
        case 4:
            return std::pow(n*pi, 4)*(std::sin(n*pi*x) + std::pow((n + 1.)/n, 4)*an*std::sin((n + 1)*pi*x) + std::pow((n + 2.)/n, 4)*bn*std::sin((n + 2)*pi*x));
        }
    }
}
