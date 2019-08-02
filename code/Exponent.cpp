double h(const std::complex<double> x, const double mu)
{
    return real(func(x, mu));
}

double H(const std::complex<double> x, const double mu)
{
    return imag(func(x, mu));
}

std::complex<double> gradient(const std::complex<double> point, const double epsilon, const double mu)
{
    return ((h(point + epsilon, mu) - h(point - epsilon, mu)) + I * (h(point + I * epsilon, mu) - h(point - I * epsilon, mu))) / (2. * epsilon);
}
