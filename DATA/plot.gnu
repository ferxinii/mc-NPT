set terminal png enhanced font 'Arial,20' size 1080,1080 enhanced

N = 30
N_atoms = 512

# LOAD P_REF INTO ARRAY
array P_ref[N]
stats 'results/pressure_ref.dat' using (P_ref[int($1)] = $2) nooutput



# PRESSURE
set output "results/pressure.png"
set title "{/:Bold Reference pressure VS MC mean virial pressure}" font ", 24"
set ylabel "P*"
set xlabel "P\\\_ref*"
set grid
set key left top box spacing 2

array P_m[N]
array P_std[N]
do for [i=1:N] {
  stats 'results/pressure_'.i.'.dat' using 5 name 'A' nooutput
  P_m[i] = A_mean
  P_std[i] = A_ssd
}
plot sample [i=1:N] '+' u (P_ref[i]):(P_m[i]):(P_std[i]) w errorbars title "MC virial P" ps 2 pt 7 lw 4 lc 6, \
     [i=1:N] '+' using (P_ref[i]):(P_ref[i]) w points title "Reference P" ps 4 lw 3 pt 2 lc 7 



# ENERGY VS PRESSURE
set output "results/energy.png"
set title "{/:Bold Energy VS reference pressure}" font ", 24"
set ylabel "<U*>"
set xlabel "P\\\_ref*"
set grid
unset key

array U_m[N]
array U_std[N]
do for [i=1:N] {
  stats 'results/energy_'.i.'.dat' using 2 name 'A' nooutput
  U_m[i] = A_mean
  U_std[i] = A_ssd
}

plot sample [i=1:N] '+' u (P_ref[i]):(U_m[i]):(U_std[i]) w errorbars notitle pt 0 lw 4 lc 6, \
     [i=1:N] '+' using (P_ref[i]):(U_m[i]) w points notitle ps 3 lw 1 pt 2 lc 2



# BOXLENGTH VS PRESSURE
set output "results/boxlength.png"
set title "{/:Bold Boxlength VS reference pressure}" font ", 24"
set ylabel "<L*>"
set xlabel "P\\\_ref*"
set yrange [6:20]
set grid
unset key

array L_m[N]
array L_std[N]
do for [i=1:N] {
  stats 'results/boxlength_'.i.'.dat' using 2 name 'A' nooutput
  L_m[i] = A_mean
  L_std[i] = A_ssd
}

plot sample [i=1:N] '+' u (P_ref[i]):(L_m[i]):(L_std[i]) w errorbars notitle pt 0 lw 4 lc 6, \
     [i=1:N] '+' using (P_ref[i]):(L_m[i]) w points notitle ps 3 lw 1 pt 2 lc 2



# DENSITY VS PRESSURE
set output "results/density.png"
set title "{/:Bold Density VS reference pressure}" font ", 24"
set ylabel "<rho*>"
set xlabel "P\\\_ref*"
set yrange [0:1]
set grid
unset key

array rho_m[N]
array rho_std[N]
do for [i=1:N] {
  stats "< LC_NUMERIC=C awk 'NR>1 {print ".N_atoms."/($2)^3}' results/boxlength_".i.".dat" name 'A' nooutput
  rho_m[i] = A_mean
  rho_std[i] = A_ssd
}

plot sample [i=1:N] '+' u (P_ref[i]):(rho_m[i]):(rho_std[i]) w errorbars notitle pt 0 lw 4 lc 6, \
     [i=1:N] '+' using (P_ref[i]):(rho_m[i]) w points notitle ps 3 lw 1 pt 2 lc 2



# ENERGY VS DENSITY
set output "results/e_vs_density.png"
set title "{/:Bold Energy VS density}" font ", 24"
set ylabel "<U*>"
set xlabel "<rho*>"
unset yrange
set grid
unset key

plot sample [i=1:N] '+' u (rho_m[i]):(U_m[i]):(U_std[i]) w errorbars notitle pt 0 lw 4 lc 6, \
     [i=1:N] '+' using (rho_m[i]):(U_m[i]) w points notitle ps 3 lw 1 pt 2 lc 2



# PRESSURE VS DENSITY (EQUATION OF STATE)
set output "results/p_vs_density.png"
set title "{/:Bold Pressure VS density}" font ", 24"
set ylabel "<P*>"
set xlabel "<rho*>"
unset yrange
set grid
unset key

plot sample [i=1:N] '+' u (rho_m[i]):(P_m[i]):(P_std[i]) w errorbars notitle pt 0 lw 4 lc 6, \
     [i=1:N] '+' using (rho_m[i]):(P_m[i]) w points notitle ps 3 lw 1 pt 2 lc 2



# GDR (RADIAL DISTRIBUTION FUNCTION)
set output "results/gdr.png"
set title "{/:Bold Radial distribution function for various pressures}" font ", 24"
set ylabel "g(r*)"
set xlabel "r*"
unset yrange
set grid
set key right top box spacing 1.5 font ", 17" maxrows 15
unset colorbox

set palette defined (0 '#3498eb', 1 '#ff3300')
plot for [i=1:N] "< awk '(NR>1){print;}' results/gdr_".i.".dat" title sprintf("P* = %.2f", P_ref[i]) with lines lc palette frac (i-1.0)/(N-1) lw 3

