set format x "% h" 
set format y "10^{%T}" 
set xlabel "Eb/No (dB)" 
set ylabel "Taxa de Erro de Bit" 
set grid
set yrange [ 1e-05 : 1.00000 ] 
set logscale y
plot  'n1024k512.dat' w linesp  t "N=1024 (AWGN)", 'saida' w l
pause -1


set terminal jpeg 12 lw 1.2 dl 1.2
set output 'sc_awgn.jpg' 
replot

set terminal postscript portrait enhanced 20 color
set size 1.6,0.7
set size ratio 0.7
set output "sc_awgn.eps"
replot

