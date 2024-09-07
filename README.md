# mc-NPT
Monte Carlo simulations of an Isobaric-Isothermal Lennard-Jones system (constant N, P, T) in Fortran.

## Physics
The partition function of the ensemble is:
$$ Q(N,P,T) = \frac{\beta P}{\lambda ^{3N}N!}\int dV V^N e^{-\beta PV} \int {d\vec{s}}^N e^{-\beta U(\vec{s}^N;L)} = \beta P\int dV e^{-\beta PV} Q (N,V,T) $$

- The code uses the Metropolis method to sample the configuration space of the system. For more details regarding reduced units and the computation of relevant quantities such as the energy and pressure, see *questions.pdf*.

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


