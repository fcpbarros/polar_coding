set format x "% h" 
set format y "10^{%T}" 
set xlabel "Eb/No (dB)" 
set ylabel "Bit Error Rate" 
set grid
set yrange [ 1e-05 : 1.00000 ] 
set logscale y
plot  'n128k64.dat' w linesp  t "N=128 and K=64",\
     'n256k128.dat' w linesp  t "N=256 and K=128",\
     'n512k256.dat' w linesp  t "N=512 and K=256",\
	  'n1024k512.dat' w linesp  t "N=1024 and K=512", \
	  'n2048k1024.dat' w linesp  t "N=2048 and K=1024"
pause -1
set terminal postscript portrait enhanced 22 color linewidth 2
set size 1.6,0.7
set size ratio 0.7
set output "sc_rate0_5.eps"
replot

