# mc-NPT
Monte Carlo simulations of an Isobaric-Isothermal Lennard-Jones system (constant N, P, T) in Fortran.

## Physics
The partition function of the ensemble is:
```math
Q(N,P,T) = \frac{\beta P}{\lambda ^{3N}N!}\int dV V^N e^{-\beta PV} \int {d\vec{s}}^N e^{-\beta U(\vec{s}^N;L)} = \beta P\int dV e^{-\beta PV} Q (N,V,T)
```

The code uses the Metropolis method to sample the configuration space of the system. For more details regarding Monte Carlo simulations, the Fortran implementation, reduced units and the computation of relevant quantities such as the energy and pressure, see *questions.pdf*.

## Compile and run
```console
[ferxinii@mb ~/mc-NPT]$ make
[ferxinii@mb ~/mc-NPT]$ ./main
```
The executable *mc_sampling* is called by *main*, and is responsible for sampling a single configuration of (N, P, T) by reading the configuration file *config.dat*. The executable *main* is responsible for managing the sweep of parameters explored (Varying P), updating *config.dat*, calling *mc_sampling* and saving the correct results.

Results are saved in *results/*. Auxiliary files are stored in *tmp/*.

To produce the plots, once *main* has run:
```console
[ferxinii@mb ~/mc-NPT]$ gnuplot plot.gnu
```

## Some results
The following are the results for the simulation of a system of Argon particles ()

<p align="center">
  <img src="/example/pressure.png" alt="Pressure" style="width: 400px;">
  <img src="/example/gdr.png" alt="Radial distribution function" style="width: 400px;">
</p>
<p align="center">
  <img src="/example/e_vs_density.png" alt="E vs rho" style="width: 400px;">
  <img src="/example/p_vs_density.png" alt="P vs rho" style="width: 400px;">
</p>
<p align="center">
  <img src="/example/pressure.png" alt="Pressure" style="width: 400px;">
  <img src="/example/energy.png" alt="Energy" style="width: 400px;">
</p>
