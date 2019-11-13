void writeB(std::vector<simplex> simplices, std::string fileName)
{
    std::ofstream file; file.open(fileName, std::ios::binary);
    if(file.is_open())
    {
        for(int index = 0; index < simplices.size(); index++)
        {
            std::complex<double> p0 = simplices[index].p0();
            std::complex<double> p1 = simplices[index].p1();
            
            if(simplices[index].isActive())
            {
                double g1 = real(p0);
                double g2 = imag(p0);
                double g3 = real(p1);
                double g4 = imag(p1);
                
                file.write((char*) &g1, sizeof(double));
                file.write((char*) &g2, sizeof(double));
                file.write((char*) &g3, sizeof(double));
                file.write((char*) &g4, sizeof(double));
            }
        }
        file.close();
    } else
    {
        std::cout << "Could not open " << fileName << std::endl;
    }
    
}

void writeB(std::vector<std::complex<double> > &psi, std::string fileName)
{
    std::ofstream file; file.open(fileName, std::ios::binary);
    if(file.is_open())
    {
        for(int index = 0; index < psi.size(); index++)
        {
            double g1 = real(psi[index]);
            double g2 = imag(psi[index]);
            
            file.write((char*) &g1, sizeof(double));
            file.write((char*) &g2, sizeof(double));
        }
    } else
    {
        std::cout << "Could not open " << fileName << std::endl;
    }
}
