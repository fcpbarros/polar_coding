set format x "% h" 
set format y "10^{%T}" 
set xlabel "Eb/No (dB)" 
set ylabel "Taxa de Erro de Bit" 
set grid
set yrange [ 1e-05 : 1.00000 ] 
set logscale y
plot 'saida2' w linesp  t "N=8 and K=4 (Rayleigh-Francisco)",\
     'ray_n8k4.dat' w linesp  t "N=8 and K=4 (Rayleigh)",\
	  'saida' w linesp  t "N=8 and K=4 (AWGN-Francisco)",\
	  'n8k4.dat' w linesp  t "N=8 and K=4 (AWGN)"
pause -1


set terminal jpeg 12 lw 1.2 dl 1.2
set output 'sc_frozen.jpg' 
replot

set terminal postscript portrait enhanced 20 color
set size 1.6,0.7
set size ratio 0.7
set output "sc_frozen.eps"
replot

