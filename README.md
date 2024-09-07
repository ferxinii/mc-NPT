# mc-NPT
Monte Carlo simulations of an Isobaric-Isothermal Lennard-Jones system (constant N, P, T) in Fortran.

## Compile and run
```console
ferxinii:~/mc-NPT$ make
ferxinii:~/mc-NPT$ ./main
```
The executable *mc_sampling* is called by *main*, and is responsible for sampling a single configuration of (N, P, T). Executable *main* is responsible for managing the sweep of parameters explored (Varying P).

## Results and plotting
Results are saved in *results/*. Auxiliary files are stored in *tmp/*.

To produce the plots, once *main* has run:
```console
ferxinii:~/mc-NPT$ gnuplot plot.gnu
```


