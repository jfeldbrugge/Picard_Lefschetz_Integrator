#include "Header.h"

std::complex<double> func(const std::complex<double> x, const double mu, const double nu)
{
    return I * nu * (pow(x - mu, 2.) / 2. + phi(x));
}

std::complex<double> phi(const std::complex<double> x)
{
    // Full integral
    if(1)
    {
//        return 1. / (1. + pow(x, 6.));
        return 2. / (1. + pow(x, 2.));
    }
    
    // Gaussian interpolation
    if(0)
    {
        const double eps = 0.42;
        const std::vector<double> coeffs = {-0.000882586, 0.00296756, -0.00532955, 0.00919606, -0.0132397, 0.0206785, -0.0280891, 0.0440207, -0.0577544, 0.0970763, -0.121122, 0.24974, -0.267257, 0.86605, -0.0748684, 0.676696, 0.0170483, 0.676696, -0.0748684, 0.86605, -0.267257, 0.24974, -0.121122, 0.0970763, -0.0577544, 0.0440207, -0.0280891, 0.0206785, -0.0132397, 0.00919606, -0.00532955, 0.00296756, -0.000882586};
        const std::vector<double> points = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5, -2.25, -2., -1.75, -1.5, -1.25, -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4.};
        
        std::complex<double> sum = 0;
        for(int i = 0; i < points.size(); i++)
        {
            sum = sum + coeffs[i] * exp(-pow((x - points[i]) / eps, 2.));
        }
        return sum;
    }
    
    // Lorentzian interpolation
    if(0)
    {
        const double eps = 0.76;
        const std::vector<double> coeffs = {-0.0250235, 0.0200499, -0.0243758, 0.00946094, -0.0203005, 0.00754246, -0.0317951, 0.0217359, -0.0749407, 0.0828628, -0.260861, 0.573068, -2.31691, 3.46795, -1.88581, 1.37217, -0.86787, 1.37217, -1.88581, 3.46795, -2.31691, 0.573068, -0.260861, 0.0828628, -0.0749407, 0.0217359, -0.0317951, 0.00754246, -0.0203005, 0.00946094, -0.0243758, 0.0200499, -0.0250235};
        const std::vector<double> points = {-4., -3.75, -3.5, -3.25, -3., -2.75, -2.5, -2.25, -2., -1.75, -1.5, -1.25, -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4.};
        std::complex<double> sum = 0;
        for(int i = 0; i < points.size(); i++)
        {
            sum = sum + coeffs[i] / (1. + pow((x - points[i]) / eps, 2.));
        }
        return sum;
    }
    
    if(0)
    {
        const double L = 2.;
        const std::vector<double> coefSin = {0, 0, 0, 0 ,0};
        const std::vector<double> coefCos = {0.1, 0, 0, 0 ,0};
        
        std::complex<double> sum = 0;
        for(int i = 0; i < min(coefSin.size(), coefCos.size()); i++)
        {
            sum = sum + coefSin[i] * sin(x * 2. * pi * double(i + 1) / L) +  coefCos[i] * cos(x * 2. * pi * double(i + 1) / L);
        }
        return sum;
    }
}

int main(int argc, const char * argv[])
{
    std::cout << "Start Integral" << std::endl;
    double nu = 10;
    
    if( argc >= 2)
    {
        nu = atof(argv[1]);
    }
    
    std::cout << "nu = " << nu << std::endl;
    
    // Parameters
    const double xMin = -20.;
    const double xMax = +20.;
    const double delta = 0.25;
    
    const double step = 0.05;
    const double thres = -20.;
    
    const int Niterations = 50;
    
    const double muMin = -4.;
    const double muMax = +4.;
    const int NMu = pow(2, 6);
    
    const int N = 6;
    const int NM = pow(2, 10);
    
    const std::string directory = "files/";
    
    // Flow of the thimble
    if(0)
    {
        const double mu = 0;
        
        std::cout << "mu = " << mu << std::endl;
        
        std::vector<simplex> simplices;
        std::vector<cp> points;
        
        initialize(simplices, points, xMin, xMax, delta);
    
        writeB(simplices, points, directory + "simplices" + std::to_string(0) + ".bin");
        for(int i = 0; i < Niterations; i++)
        {
            flow(points, step, mu, nu, thres);
            subdevide(simplices, points, delta);
            clean(simplices);
            writeB(simplices, points, "files/simplices" + std::to_string(i + 1) + ".bin");
            std::cout << i << " Integral = " << integrate(simplices, points, mu, nu, N) << std::endl;
        }
    }
    
    // Compute the Picard-Lefschetz thimble
    if(1)
    {
        std::vector<double> mus(NMu);
        externalVariables(mus, muMin, muMax, NMu);
        
        std::cout << "Picard-Lefschetz:";
        std::vector<std::vector<simplex>> PL(mus.size());
#pragma omp parallel for
        for(int index = 0; index < mus.size(); index++)
        {
            std::cout << "."; std::cout.flush();
            std::vector<simplex> simplices;
            std::vector<cp> points;
            initialize(simplices, points, xMin, xMax, delta);
            flow(simplices, points, step, mus[index], nu, thres, delta, Niterations);
            clean(simplices);
            importPoints(simplices, points);
            PL[index] = simplices;
        }
        std::cout << std::endl;
        
        writeB(PL, mus,
               "files/PL_nu=" + std::to_string(int(nu)) + ".bin",
               "files/log_nu=" + std::to_string(int(nu)) + ".txt",
               "files/mus_nu=" + std::to_string(int(nu)) + ".bin");
        
        
        // Evaluate the integral
        std::cout << "Integrate:";
        std::vector<std::complex<double>> result(NM);
        
#pragma omp parallel for
        for(int iMu = 0; iMu < NM; iMu++)
        {
            if((iMu % (NM / 16)) == 0){std::cout << "."; std::cout.flush();}
            const double mu = (muMax - muMin) / (NM - 1) * iMu + muMin;
            int index = nearest(mus, mu);
            result[iMu] = integrate(PL[index], mu, nu, N);
        }
        std::cout << std::endl;
        
        writeB(result, "files/result_nu=" + std::to_string(int(nu)) + ".bin");
    }
    
    std::cout << "Done" << std::endl << std::endl;
    return 0;
}
