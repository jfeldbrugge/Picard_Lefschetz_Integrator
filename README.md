# Picard-Lefschetz Integrator

Oscillatory integrals, such as the ones occurring in optical interference phenomena, are central to physics. These integrals are however often expensive to evaluate numerically. Picard-Lefschetz theory can be used to transform analytic oscillatory integrals into sums of convex integrals by deforming the integration domain in the complex plane. In the paper ..., we propose a new numerical integrator, based on Picard-Lefschetz theory, which enables us the efficiently evaluate the interference effects near caustics in lenses. We here present an accompanying code for the one-dimensional integral.

## Picard-Lefschetz theory ##
We start with the original integration domain and continuously deform the domain in the complex plane using the downward flow of the exponent. Consider for example the integral

<img src="figures/Integral.png" height="75" />

with the exponent

<img src="figures/Exponent.png" height="75" />

with the strength of the lense <img src="figures/alpha.png" height="15" /> frequency <img src="figures/nu.png" height="15" />, as a function of <img src="figures/mu.png" height="15" />. We write the integral as a sum of integrals along the steepest descent contours

<img src="figures/Integral2.png" height="75" />

The trick is to find an approximation of the steepest descent contours <img src="figures/Descent.png" height="15" /> before evaluating the integral.

We first expand the exponent into a real and an imaginary part

<img src="figures/RealImag.png" height="40" />

The real part is known as the Morse function. The downward flow is defined as 

<img src="figures/FlowEquation.png" height="75" />

with the boundary condition

<img src="figures/BoundaryCondition.png" height="40" />

The flow induces a deformation of the original integration contour

<img src="figures/Deformation.png" height="40" />

which, in the limit of large <img src="figures/lambda.png" height="15" />, approaches the thimble 

<img src="figures/Limit.png" height="40" />

The deformation is illustrated in the figure:

<img src="figures/Flow.gif" height="600" />

where the red points are the saddle points, the black curves are the steepest descent and ascent curves and the blue line is the deformed integration domain. Along the thimble, the evaluation of the integral is trivial as the integrand is convex. This leads to an interference pattern, in which we can observe two single image region and a triple image region, separated by a fold caustic. See the figure for the absolute square of the wavefunction, representing the intensity of the lensed image

<img src="figures/Interference.png" height="300" />

For more details see the paper ...

## Implementation ##
The code is written in C++ using only the standard libraries. Compile the code with

```console
$ g++ main.cpp -o PicardLefschetz.a
```

Evaluate the code by running the executable 

```console
$ ./PicardLefschetz.a
```

The program outputs the binary files "simplices#.bin", "psi.bin". The files "simplices#.bin" are thimbles. The file "Psi.bin" is the integral as a function of <img src="figures/mu.png" height="15" />.

## Collaborators ##

* Job Feldbrugge
* Neil Turok 
* Ue-Li Pen
