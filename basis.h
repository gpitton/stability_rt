#ifndef BASIS_HEADER
#define BASIS_HEADER


constexpr double pi = 3.141592653589793238463;


enum parity { Even, Odd };


class BasisFunction {
 public:
    virtual double operator()(const double) = 0;
    virtual double deriv(int, const double) = 0;
};


class Trig3 : public BasisFunction {
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
