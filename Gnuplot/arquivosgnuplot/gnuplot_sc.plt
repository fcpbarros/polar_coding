set format x "% h" 
set format y "10^{%T}" 
set xlabel "Eb/No (dB)" 
set ylabel "Taxa de Erro de Bit" 
set grid
set yrange [ 1e-05 : 1.00000 ] 
set logscale y
plot 'ray_n128k64.dat' w linesp  t "N=128 (Rayleigh)",\
     'ray_n256k128.dat' w linesp  t "N=256 (Rayleigh)",\
     'ray_n512k256.dat' w linesp  t "N=512 (Rayleigh)",\
     'ray_n1024k512.dat' w linesp  t "N=1024 (Rayleigh)",\
     'ray_n2048k1024.dat' w linesp  t "N=2048 (Rayleigh)",\
     'n128k64.dat' w linesp  t "N=128 (AWGN)",\
     'n256k128.dat' w linesp  t "N=256 (AWGN)",\
     'n512k256.dat' w linesp  t "N=512 (AWGN)",\
	  'n1024k512.dat' w linesp  t "N=1024 (AWGN)", \
	  'n2048k1024.dat' w linesp  t "N=2048 (AWGN)"
pause -1


set terminal jpeg 12 lw 1.2 dl 1.2
set output 'sc_rate0_5.jpg' 
replot

set terminal postscript portrait enhanced 20 color
set size 1.6,0.7
set size ratio 0.7
set output "sc_rate0_5.eps"
replot

