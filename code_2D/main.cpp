#include "Header.h"

std::complex<double> func(const pointC p, const pointD mu, const double nu)
{
    return I * nu * (pow(p.x - mu.x, 2.) / 2. + pow(p.y - mu.y, 2.) / 2. + phi(p));
}

std::complex<double> phi(const pointC p)
{
    // Full integral
    if(1)
    {
        return 0.65 / (1. + p.x * p.x + 2. * p.y * p.y);
    }

    // Gaussian interpolation
    if(0)
    {
        const double eps = 1;
        const std::vector<double> coeffs = {1};
        const std::vector<pointD> points = {pointD(0,0)};

        std::complex<double> sum = 0;
        for(int i = 0; i < points.size(); i++)
        {
            sum = sum + coeffs[i] * exp(-pow((p.x - points[i].x) / eps, 2.) - pow((p.y - points[i].y) / eps, 2.));
        }
        return sum;
    }

    // Lorentzian interpolation
    if(0)
    {
        const double eps = 1;
        const std::vector<double> coeffs = {1};
        const std::vector<pointD> points = {pointD(0,0)};
        
        std::complex<double> sum = 0;
        for(int i = 0; i < points.size(); i++)
        {
            sum = sum + coeffs[i] / (1. + pow((p.x - points[i].x) / eps, 2.) + pow((p.y - points[i].y) / eps, 2.));
        }
        return sum;
    }

    // Fourier
    if(0)
    {
        const double L = 2.;
        const std::vector<double> coefSin = {0., -0.00107114, -0.000623854, 0.00118156, -0.000196488, -0.000509684, 0.0012452, -0.000158929, 0.000632952, -0.0001544, -0.0000485848, -0.000462906, -0.00118477, 0.000640332, -0.0000335514, 0.00028905, 0.00032202, -0.000529926, -0.000938008, -0.00183425, -0.000678926, 0.000269456, -0.0000328556, -0.000623054, 0.000702522, 0.00053917, -0.000044644, -0.000582282, 0.00015913, 0.000210624, -0.000528754, -0.000228756, -0.000126658, -0.000242426, -0.000205752, -0.000024445};
        const std::vector<double> coefCos = {0., 0.00176056, -0.00113817, 0.000263744, 0.000063046, 0.000308638, 0.000268442, 0.000474586, -0.000498962, 0.000196472, -0.000254158, 0.000546198, -0.00103912, 0.000902562, -0.0000232504, 0.00137611, 0.000390128, 0.0000079428, -0.000439484, 0.00058306, 0.000542024, 0.000779034, -0.000788798, 0.0000974324, -0.000742616, 0.000198539, -0.000168164, -0.000129533, -0.000400742, 0.000359572, 0.000430736, 0.000153811, -0.0000603164, 0.000465384, 0.0000479118, -0.000275958};
        const std::vector<pointD> points = {
            pointD(0,0), pointD(0,1), pointD(0,2), pointD(0,3), pointD(0,4), pointD(0,5),
            pointD(1,0), pointD(1,1), pointD(1,2), pointD(1,3), pointD(1,4), pointD(1,5),
            pointD(2,0), pointD(2,1), pointD(2,2), pointD(2,3), pointD(2,4), pointD(2,5),
            pointD(3,0), pointD(3,1), pointD(3,2), pointD(3,3), pointD(3,4), pointD(3,5),
            pointD(4,0), pointD(4,1), pointD(4,2), pointD(4,3), pointD(4,4), pointD(4,5),
            pointD(5,0), pointD(5,1), pointD(5,2), pointD(5,3), pointD(5,4), pointD(5,5)};

        std::complex<double> sum = 0;
        for(int i = 0; i < min(coefSin.size(), coefCos.size(), points.size()); i++)
        {
            sum = sum + (coefSin[i] * sin((p.x * points[i].x + p.y * points[i].y) * 2. * pi / L) +
                         coefCos[i] * cos((p.x * points[i].x + p.y * points[i].y) * 2. * pi / L));
        }
        return 5. * sum;
    }
}

int main(int argc, const char * argv[])
{
    std::cout << "Start Integral" << std::endl;
    double nu = 50.;
    
    if( argc >= 2)
    {
        nu = atof(argv[1]);
    }
    
    std::cout << "nu = " << nu << std::endl;
    
    // Parameters
    const double xMin = -5.;
    const double xMax = +5.;
//        const double xMin = -2.;
//        const double xMax = +2.;
//    const double delta = 0.25;
    const double delta = 0.1;
//    const double delta = 0.05;
    
    const double step = 0.02;
    const double thres = -20.;
    
    const int Niterations = 20;
//    const int Niterations = 15;
    
//    const double muMin = -.5;
//    const double muMax = +.5;
//    const double muMin = -1.;
//    const double muMax = +1.;
    
    const double muMin = -4.;
    const double muMax = +4.;

    const int NMu = pow(2, 6);
//    const int NMu = pow(2, 4);
    
    const int N = 3;
//    const int N = 4;
//    const int NM = pow(2, 7);
//    const int NM = pow(2, 8);
    const int NM = pow(2, 12);
    const std::string directory = "files/";
    
    // Flow of the thimble
    if(0)
    {
        const pointD mu(0, 0);
        std::cout << "mu = (" << mu << ")" << std::endl;
        
        std::vector<simplex> simplices;
        std::vector<cp> points;
        
        initialize(simplices, points, xMin, xMax, delta);
        
        writeB(simplices, points, directory + "simplices" + std::to_string(0) + ".bin");
        for(int i = 0; i < Niterations; i++)
        {
            flow(points, step, mu, nu, thres);
            subdevide(simplices, points, delta);
            writeB(simplices, points, "files/simplices" + std::to_string(i + 1) + ".bin");
            std::cout << i << " Integral = " << integrate(simplices, points, mu, nu, N) << std::endl;
        }
    }
    
    // Compute the Picard-Lefschetz thimble
    if(1)
    {
        std::vector<pointD> mus(NMu * NMu);
        externalVariables(mus, muMin, muMax, NMu);
        
        std::cout << "Picard-Lefschetz:";
        std::vector<std::vector<simplex>> PL(mus.size());
#pragma omp parallel for
        for(int index = 0; index < mus.size(); index++)
        {
            if((index % (NMu * NMu / 16)) == 0){std::cout << "."; std::cout.flush();}
            std::vector<simplex> simplices;
            std::vector<cp> points;
            initialize(simplices, points, xMin, xMax, delta);
            flow(simplices, points, step, mus[index], nu, thres, delta, Niterations);
            importPoints(simplices, points);
            PL[index] = simplices;
        }
        std::cout << std::endl;
        
//        writeB(PL, mus,
//               "files/PL_nu=" + std::to_string(int(nu)) + ".bin",
//               "files/log_nu=" + std::to_string(int(nu)) + ".txt",
//               "files/mus_nu=" + std::to_string(int(nu)) + ".bin");

        // Evaluate the integral
        std::cout << "Integrate:       ";
        std::vector<std::complex<double>> result(NM * NM);

#pragma omp parallel for
        for(int iMu = 0; iMu < NM; iMu++)
        {
            if((iMu % (NM / 16)) == 0){std::cout << "."; std::cout.flush();}
            
            for(int jMu = 0; jMu < NM; jMu++)
            {
                const pointD mu((muMax - muMin) / (NM - 1) * iMu + muMin,
                                (muMax - muMin) / (NM - 1) * jMu + muMin);
                int index = nearest(mus, mu);
                result[jMu + iMu * NM] = integrate(PL[index], mu, nu, N);
            }
        }
        std::cout << std::endl;

        writeB(result, "files/result_nu=" + std::to_string(int(nu)) + ".bin");
    }
    
    std::cout << "Done" << std::endl << std::endl;
    return 0;
}
