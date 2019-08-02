#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <fstream>

const double epsilon = std::numeric_limits<double>::epsilon() ;
const double pi = 3.141592654;
const std::complex<double> I (0,1.0);

class simplex {
private:
    std::complex<double> pnts[2];
    bool active;
    
public:
    simplex(std::complex<double> p0, std::complex<double> p1);
    void flow(const double tau, const double mu, const double xMin, const double xMax, const double thres);
    void subdivide(std::vector<simplex> &simplices, const double delta);
    void refine(std::vector<simplex> &simplices, const double mu, const double nu, const double accuracy);
    std::complex<double> integrate(const double mu, const double nu);
    bool isActive();
    std::complex<double> p0();
    std::complex<double> p1();
};

std::complex<double> func(const std::complex<double> x, const double mu);
double h(const std::complex<double> x, const double mu);
double H(const std::complex<double> x, const double mu);
std::complex<double> gradient(const std::complex<double> point, const double epsilon, const double mu);

void initialize(std::vector<simplex> &simplices, const double xMin, const double xMax, const double delta);
void clean(std::vector<simplex> &simplices);
void subdivide(std::vector<simplex> &simplices, const int index, const double delta);
void refine(const std::vector<std::vector<simplex> > &PL, std::vector<std::vector<simplex> > &PLrefined, const double muMin, const double deltaMu, const double nu, const double accuracy);
void flow(std::vector<simplex> &simplices, const int Niterations, const double tau, const double mu, const double xMin, const double xMax, const double delta, const double thres);

std::complex<double> integrate(std::vector<simplex> &simplices, const double mu, const double nu);

void writeB(std::vector<simplex> simplices, std::string fileName);
void writeB(std::vector<std::complex<double> > &psi, std::string fileName);

#include "Exponent.cpp"
#include "Simplex.cpp"
#include "Utility.cpp"
