# Picard-Lefschetz Integrator
Oscillatory integrals, such as the ones occurring in optical interference phenomena, are central to physics. These integrals are however often expensive to evaluate numerically. Picard-Lefschetz theory can be used to transform analytic oscillatory integrals into sums of convex integrals by deforming the integration domain in the complex plane. In the paper <a href="url">Feldbrugge, Pen and Turok, 2019</a>, we propose a new numerical integrator, which enables us the efficiently evaluate the interference effects near caustics in lenses. We here present corresponding code for the one-dimensional integral. The code serves as a proof of concept and is not optized. Different applications might require different implementations.

## Picard-Lefschetz theory ##
We illustrate the technique by studying the one-dimensional example

<img src="figures/Integral.png" height="75" />

with the exponent

<img src="figures/Exponent.png" height="75" />

with the strength of the lense <img src="figures/alpha.png" height="15" />, the frequency <img src="figures/nu.png" height="15" />, as a function of <img src="figures/mu.png" height="15" />. 

Using Picard-Lefschetz theory we can express the integral as a sum of convex integrals along the steepest descent contours <img src="figures/Descent.png" height="15" /> of the exponent, i.e.,

<img src="figures/Integral2.png" height="75" />

The evaluation of the integral is thus reduced to the finding the relevant steepest descent contours. We first expand the exponent into a real and an imaginary part

<img src="figures/RealImag.png" height="40" />

The downward flow, corresponding to the real part of the exponent, is defined as 

<img src="figures/FlowEquation.png" height="75" />

with the boundary condition

<img src="figures/BoundaryCondition.png" height="40" />

The flow induces a deformation of the original integration contour

<img src="figures/Deformation.png" height="40" />

which, in the limit of large <img src="figures/lambda.png" height="15" />, approaches the thimble 

<img src="figures/Limit.png" height="40" />

The deformation for $\mu=0$ is illustrated in the figure

<img src="figures/Flow.gif" height="600" />

where the red points are the saddle points, the black curves are the steepest descent and ascent contours eminated by the saddle points, the blue curve is the deformed integration domain consisting of a number of line segments and the grey areas are the regions for which the real part of the exponent is smaller than a threshold. Along the thimble, the evaluation of the integral is trivial as the integrand is convex. 

The resulting interference pattern consists of a single image region in which the intensity is small, a triple-image region in which the intensity oscillates and two fold caustics separating the single- and triple-image regions.

<img src="figures/Interference.png" height="300" />

For more details see the paper <a href="url">Feldbrugge, Pen and Turok, 2019</a> and the wiki.

## Implementation ##
The code is written in C++. The implementation relies only the standard libraries. Compile the code with

```console
$ g++ main.cpp -o PicardLefschetz.a
```

Evaluate the integral described above by running the executable 

```console
$ ./PicardLefschetz.a
```

The program outputs the binary files "simplices#.bin" consisting the of the line-segments in the thimble, and "psi#.bin" consisting of the integral for various <img src="figures/mu.png" height="15" />. For more information see the wiki.

## Collaborators ##

* Job Feldbrugge
* Neil Turok 
* Ue-Li Pen
