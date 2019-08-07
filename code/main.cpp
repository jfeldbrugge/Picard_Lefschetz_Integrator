#include "Header.h"

std::complex<double> func(const std::complex<double> x, const double mu)
{
    const double alpha = 2.;
    return I * (pow(x - mu, 2) + alpha / (1. + x * x));
}

int main(int argc, const char * argv[]) {
    std::cout << "Start Integral" << std::endl;
    
    // Parameters
    const double xMin = -4.;
    const double xMax = +4.;
    const double delta = 0.25;

    const double tau = 0.01;
    const double thres = -2.;
    const double accuracy = 0.05;
    
    const double muMin = -0.4;
    const double muMax = 0.4;
    const double deltaMu = 0.05;
    std::vector<double> nu_list = {50, 100, 500};

    // Flow of the thimble
    if(1)
    {
        const double mu = 0.;
        std::vector<simplex> simplices;
        initialize(simplices, xMin, xMax, delta);
        writeB(simplices, "simplices" + std::to_string(0) + ".bin");
        
        for(int i = 0; i < 150; i++)
        {
            for(int index = 0; index < simplices.size(); index++)
            {
                simplices[index].flow(tau, mu, xMin, xMax, thres);
                simplices[index].subdivide(simplices, delta);
            }
            
            clean(simplices);
            writeB(simplices, "simplices" + std::to_string(i + 1) + ".bin");
        }
        writeB(simplices, "simplices.bin");
    }
    
    // Perform integral for various mu and nu
    if(1)
    {
        std::vector<std::vector<simplex> > PL;

        // Evaluate course thimble
        for(double mu = muMin; mu < muMax + epsilon; mu = mu + deltaMu)
        {
            std::vector<simplex> simplices;
            initialize(simplices, xMin, xMax, delta);
            
            // Flow the integration domain
            std::cout << mu << " ";
            flow(simplices, 150, tau, mu, xMin, xMax, delta, thres);
            PL.push_back(simplices);
            
            // Output flowed integration domain
            writeB(simplices, "simplices_mu=" + std::to_string(round(mu * 100) / 100).substr(0,5) + "_nu=" + std::to_string(1) + ".bin");
        }
        std::cout << std::endl << std::endl;
        
        // Iterate over the frequencies
        for(int nuElement = 0; nuElement < nu_list.size(); nuElement++)
        {
            const double nu = nu_list[nuElement];
            
            // Refine thimble for this nu
            std::cout << "nu = " << nu << std::endl;
            std::vector<std::vector<simplex> > PLrefined;
            refine(PL, PLrefined, muMin, deltaMu, nu, accuracy);
            
            // Integrate for this nu over refined thimble
            std::vector<std::complex<double> > psi;
            for(double mu = muMin; mu < muMax + 0.00000001; mu = mu + 0.0001)
            {
                const int index = int((mu - muMin) / deltaMu + 0.000001);
                psi.push_back(integrate(PLrefined[index], mu, nu));
            }
            
            // Output refined thimble
            for(int index = 0; index < PLrefined.size(); index++)
            {
                const double mu = index * deltaMu + muMin;
                std::vector<simplex> simplices = PLrefined[index];
                writeB(simplices, "simplices_mu=" + std::to_string(round(mu * 100) / 100).substr(0,5) + "_nu=" + std::to_string(int(nu)) + ".bin");
            }

            // Output integral over refined thimble 
            writeB(psi, "psi_nu=" + std::to_string(int(nu)) + ".bin");
        }
    }

    std::cout << "Done" << std::endl;
}
