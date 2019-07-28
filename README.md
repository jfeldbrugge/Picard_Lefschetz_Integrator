# Picard-Lefschetz Integrator
Oscillatory integrals, such as the ones occurring in optical interference phenomena, are central to physics. These integrals are however often expensive to evaluate numerically. Picard-Lefschetz theory can be used to transform analytic oscillatory integrals into sums of convex integrals by deforming the integration domain in the complex plane. In the paper <a href="url">Feldbrugge, Pen, and Turok 2019</a>, we propose a new numerical integrator, which enables us the efficiently evaluate the interference effects near caustics in lenses. We here present corresponding code for the one-dimensional integral. The code serves as a proof of concept and is not optimal. Different applications might require different implementations.

<img src="figures/2DCaustics.png" height="400" />

## Picard-Lefschetz theory ##
We illustrate the technique by studying the one-dimensional example

<a href="https://www.codecogs.com/eqnedit.php?latex=\Psi(\mu;\nu)=\sqrt{\frac{\nu}{\pi}}\int&space;e^{i\phi(x;\mu)\nu}\mathrm{d}x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Psi(\mu;\nu)=\sqrt{\frac{\nu}{\pi}}\int&space;e^{i\phi(x;\mu)\nu}\mathrm{d}x" title="\Psi(\mu;\nu)=\sqrt{\frac{\nu}{\pi}}\int e^{i\phi(x;\mu)\nu}\mathrm{d}x" /></a>

with the exponent

<a href="https://www.codecogs.com/eqnedit.php?latex=i\phi(x;\mu)\nu&space;=&space;i&space;\left[(x-\mu)^2&space;&plus;&space;\frac{\alpha}{1&plus;x^2}&space;\right&space;]\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?i\phi(x;\mu)\nu&space;=&space;i&space;\left[(x-\mu)^2&space;&plus;&space;\frac{\alpha}{1&plus;x^2}&space;\right&space;]\nu" title="i\phi(x;\mu)\nu = i \left[(x-\mu)^2 + \frac{\alpha}{1+x^2} \right ]\nu" /></a>

with the strength of the lense <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha" title="\alpha" /></a>, the frequency <a href="https://www.codecogs.com/eqnedit.php?latex=\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu" title="\nu" /></a>, as a function of <a href="https://www.codecogs.com/eqnedit.php?latex=\mu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu" title="\mu" /></a>. 

Using Picard-Lefschetz theory we can express the integral as a sum of convex integrals along the steepest descent contours <a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{J}_i\subset&space;\mathbb{C}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{J}_i\subset&space;\mathbb{C}" title="\mathcal{J}_i\subset \mathbb{C}" /></a> of the real part of the exponent, <i>i.e.</i>,

<a href="https://www.codecogs.com/eqnedit.php?latex=\Psi(\mu;\nu)&space;=&space;\sqrt{\frac{\nu}{\pi}}&space;\sum_i&space;\int_{\mathcal{J}_i}e^{i\phi(x;\mu)\nu}\mathrm{d}x" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Psi(\mu;\nu)&space;=&space;\sqrt{\frac{\nu}{\pi}}&space;\sum_i&space;\int_{\mathcal{J}_i}e^{i\phi(x;\mu)\nu}\mathrm{d}x" title="\Psi(\mu;\nu) = \sqrt{\frac{\nu}{\pi}} \sum_i \int_{\mathcal{J}_i}e^{i\phi(x;\mu)\nu}\mathrm{d}x" /></a>

The evaluation of the integral is thus reduced to finding the relevant steepest descent contours. We first expand the exponent into a real and an imaginary part

<a href="https://www.codecogs.com/eqnedit.php?latex=i\phi(x;\mu)=h(x;\mu)&space;&plus;&space;i&space;H(x;\mu)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?i\phi(x;\mu)=h(x;\mu)&space;&plus;&space;i&space;H(x;\mu)" title="i\phi(x;\mu)=h(x;\mu) + i H(x;\mu)" /></a>

The downward flow, corresponding to the real part of the exponent <a href="https://www.codecogs.com/eqnedit.php?latex=h(x)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?h(x)" title="h(x)" /></a>, is defined as 

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\gamma_\lambda(x)}{\partial&space;\lambda}=-\nabla_x&space;h[\gamma_\lambda(x)]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\gamma_\lambda(x)}{\partial&space;\lambda}=-\nabla_x&space;h[\gamma_\lambda(x)]" title="\frac{\partial \gamma_\lambda(x)}{\partial \lambda}=-\nabla_x h[\gamma_\lambda(x)]" /></a>

with the boundary condition <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma_0(x)=x&space;\in&space;\mathbb{C}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma_0(x)=x&space;\in&space;\mathbb{C}" title="\gamma_0(x)=x \in \mathbb{C}" /></a>. The flow induces a deformation of the original integration domain <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{R}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{R}" title="\mathbb{R}" /></a>, <i>i.e.</i>,

<a href="https://www.codecogs.com/eqnedit.php?latex=X_\lambda&space;=&space;\gamma_\lambda(\mathbb{R})&space;\subset&space;\mathbb{C}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X_\lambda&space;=&space;\gamma_\lambda(\mathbb{R})&space;\subset&space;\mathbb{C}" title="X_\lambda = \gamma_\lambda(\mathbb{R}) \subset \mathbb{C}" /></a>

which, in the limit of large <a href="https://www.codecogs.com/eqnedit.php?latex=\lambda" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lambda" title="\lambda" /></a>, approaches the thimble 

<a href="https://www.codecogs.com/eqnedit.php?latex=\lim_{\lambda\to\infty}X_\lambda&space;=&space;\mathcal{J}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\lim_{\lambda\to\infty}X_\lambda&space;=&space;\mathcal{J}" title="\lim_{\lambda\to\infty}X_\lambda = \mathcal{J}" /></a>

The deformation for <a href="https://www.codecogs.com/eqnedit.php?latex=\mu=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu=0" title="\mu=0" /></a> is illustrated in the figure

<img src="figures/Flow.gif" height="600" />

The red points are the saddle points, the black curves are the steepest descent and ascent contours emanated by the saddle points, the blue curve is the deformed integration domain consisting of a number of line segments and the grey areas are the regions for which the real part of the exponent is smaller than a threshold. Along the thimble, the evaluation of the integral is trivial as the integrand is convex. 

The resulting interference patterns for various <a href="https://www.codecogs.com/eqnedit.php?latex=\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu" title="\nu" /></a> consist of a single image region in which the intensity is small, a triple-image region in which the intensity oscillates and two fold caustics separating the single- and triple-image regions, as can be observed in the figure plotting the intensity <a href="https://www.codecogs.com/eqnedit.php?latex=I&space;(\mu;\nu)=&space;|\Psi(\mu;\nu)|^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?I&space;(\mu;\nu)=&space;|\Psi(\mu;\nu)|^2" title="I (\mu;\nu)= |\Psi(\mu;\nu)|^2" /></a>.

<img src="figures/Interference.png" height="300" />

For more details on the proceedure <a href="url">Feldbrugge, Pen, and Turok 2019</a> and the [wiki](https://github.com/jfeldbrugge/Picard_Lefschetz_Integrator/wiki).

## Implementation ##
The code is written in C++. The implementation relies only on the standard libraries. Compile the code with

```console
$ g++ -std=c++11 main.cpp -o PicardLefschetz.a
```

Evaluate the oscillatory integral by running the executable 

```console
$ ./PicardLefschetz.a
```

The program outputs the to types of binary files. The thimbles are saved in the files "simplices#.bin". The integrals for various <a href="https://www.codecogs.com/eqnedit.php?latex=\mu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu" title="\mu" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\nu" title="\nu" /></a> are saved in the files "psi#.bin". For more information see the [wiki](https://github.com/jfeldbrugge/Picard_Lefschetz_Integrator/wiki).

## Collaborators ##

* Job Feldbrugge
* Neil Turok 
* Ue-Li Pen
