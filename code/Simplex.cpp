// simplex methods
simplex::simplex(std::complex<double> p0, std::complex<double> p1)
{
    pnts[0] = p0;
    pnts[1] = p1;
    active = true;
}

void simplex::flow(const double tau, const double mu,
                   const double xMin, const double xMax, const double thres)
{
    active = (xMin < real(pnts[0]) && real(pnts[0]) < xMax &&
              xMin < real(pnts[1]) && real(pnts[1]) < xMax &&
              h(pnts[0], mu) > thres && h(pnts[1], mu) > thres);
    if(active)
    {
        pnts[0] = pnts[0] - tau * gradient(pnts[0], 0.01, mu);
        pnts[1] = pnts[1] - tau * gradient(pnts[1], 0.01, mu);
    }
}

void simplex::subdivide(std::vector<simplex> &simplices, const double delta)
{
    if(active && abs(pnts[0] - pnts[1]) > delta)
    {
        active = false;
        
        simplex simp0(pnts[0], (pnts[0] + pnts[1]) / 2.);
        simplex simp1((pnts[0] + pnts[1]) / 2., pnts[1]);
        
        simplices.push_back(simp0);
        simplices.push_back(simp1);
    }
}

void simplex::refine(std::vector<simplex> &simplices, const double mu, const double nu, const double accuracy)
{
    if(active)
    {
        active = false;
        const double thres = 4. * log(accuracy);
        if(h(pnts[0], mu) * nu > thres || h(pnts[1], mu) * nu > thres)
        {
            simplex simp0(pnts[0], (pnts[0] + pnts[1] ) / 2.);
            simplex simp1((pnts[0] + pnts[1] ) / 2., pnts[1]);
            simplices.push_back(simp0);
            simplices.push_back(simp1);
        }
    }
}

std::complex<double> simplex::integrate(const double mu, const double nu)
{
    if(active)
    {
        // Trapezoidal rule
        return (std::exp(func(pnts[0], mu) * nu) + std::exp(func(pnts[1], mu) * nu)) * (pnts[1] - pnts[0]) / 2.;
        // Midpoint rule
//        return std::exp(func((pnts[0] + pnts[1]) / 2., mu) * nu) * (pnts[1] - pnts[0]);
    } else
    {
        return 0.;
    }
}

bool simplex::isActive()
{
    return active;
}

std::complex<double> simplex::p0()
{
    return pnts[0];
}

std::complex<double> simplex::p1()
{
    return pnts[1];
}

// Simplex related routines
void initialize(std::vector<simplex> &simplices, const double xMin, const double xMax, const double delta)
{
    for(double x = xMin; x < xMax + epsilon; x = x + delta)
    {
        simplex sim(x, x + delta);
        simplices.push_back(sim);
    }
}

void flow(std::vector<simplex> &simplices, const int Niterations, const double tau, const double mu, const double xMin, const double xMax, const double delta, const double thres)
{
    for(int i = 0; i < Niterations; i++)
    {
        int length = int(simplices.size());
        for(int index = 0; index < length; index++)
        {
            simplices[index].flow(tau, mu, xMin, xMax,  thres);
            simplices[index].subdivide(simplices, delta);
        }
        clean(simplices);
    }
}

void clean(std::vector<simplex> &simplices)
{
    std::vector<simplex> newSimplices;
    for(int index = 0; index < simplices.size(); index++)
    {
        simplex element = simplices[index];
        if(element.isActive())
        {
            newSimplices.push_back(element);
        }
    }
    simplices = newSimplices;
}

void refine(const std::vector<std::vector<simplex> > &PL, std::vector<std::vector<simplex> > &PLrefined, const double muMin, const double deltaMu, const double nu, const double accuracy)
{
    for(int index = 0; index < PL.size(); index++)
    {
        std::vector<simplex> simplices = PL[index];
        const double mu = index * deltaMu + muMin;
        
        std::complex<double> f1 = 2, f2 = 1;
        while (abs(f1 - f2) / abs(f2) > accuracy)
        {
            f1 = integrate(simplices, mu, nu);
            int num = int(simplices.size());
            for(int index = 0; index < num; index++)
            {
                simplices[index].refine(simplices, mu, nu, accuracy);
            }
            f2 = integrate(simplices, mu, nu);
        }
        
        clean(simplices);
        PLrefined.push_back(simplices);
    }
}

std::complex<double> integrate(std::vector<simplex> &simplices, const double mu, const double nu)
{
    std::complex<double> sum = 0;
    for(int index = 0; index < simplices.size(); index++)
    {
        sum = sum +  simplices[index].integrate(mu, nu);
    }
    
    return sqrt(nu / pi) * sum;
}
