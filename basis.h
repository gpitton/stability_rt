#ifndef BASIS_HEADER
#define BASIS_HEADER


constexpr double pi = 3.141592653589793238463;

enum parity { Even, Odd };


class Trig3 {
 public:
    Trig3(int);
    double operator()(const double);
    double deriv(int, const double);
 private:
    int      n;
    parity   t;
    double   an, bn;
};

#endif
