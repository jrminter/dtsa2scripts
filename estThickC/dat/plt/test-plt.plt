# Gnuplot
# Plots the spectrum generated by 'penepma16'
#
# 2020-01-06
# By J. R. Minter
# This version uses the genuine Helvetica.ttf font
# 
set terminal postscript enhanced color "Helvetica" 14
set output "test-plt.ps"
# proc-gnuplot test-plt
#
# set terminal png
# set terminal window
# set terminal aqua
# set terminal window

lablFont = "Helvetica,14"
titlFont = "Helvetica,16"
keyFont  = "Helvetica,14"
ticFont  = "Helvetica,12"

set key font keyFont
set key right

set tics font ticFont

unset xzeroaxis
set zero 1.0e-60
set xzeroaxis
set xrange [0:15]
set yrange [1e-11:*]

set title 'PENEPMA Simulation of EagleXG at 15 kV' font titlFont
set xlabel "energy [keV]" font lablFont
set ylabel "PDF [1/(eV*sr*electron)]" font lablFont

set logscale y

plot 'pe-spect-01.dat' u ($1/1000.):2:3 w lines lw 2 lc rgb "#0000A0" title "spectrum"
