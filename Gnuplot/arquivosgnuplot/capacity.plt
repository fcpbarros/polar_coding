set format x "% h" 
set format y "10^{%T}" 
set xlabel "Eb/No (dB)" 
set ylabel "Taxa de Erro de Bit" 
set grid
set yrange [ 1e-05 : 1.00000 ] 
set logscale y
#
SNR(x) = 10**(x/10)
#
f(x) = 0.5 * log(1+SNR(x))
plot f(x),'n128k64.dat' w linesp  t "N=128 (AWGN)",\
     'n256k128.dat' w linesp  t "N=256 (AWGN)",\
     'n512k256.dat' w linesp  t "N=512 (AWGN)",\
	  'n1024k512.dat' w linesp  t "N=1024 (AWGN)", \
	  'n2048k1024.dat' w linesp  t "N=2048 (AWGN)"
pause -1



