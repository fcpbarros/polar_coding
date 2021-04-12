function [tGauss_In, tGauss_Qn] = gaussGen(H, N, m)
%% White Gaussian noise Generation (Frequency domain)
fStdGauss_I = randn(1,N); % In-phase component
fStdGauss_Q = randn(1,N); % Quadrature component
fGauss_I = H.*fStdGauss_I; % In-phase component shaped by the ...Doppler filter
fGauss_Q = H.*fStdGauss_Q; % Quadrature component shaped by the ...Doppler filter

%% Evaluating Inverse Fast Fourier Transform (Time domain)
tGauss_I = real(ifft(fGauss_I,N)); % Correlated in-phase Gaussian ...vector
tGauss_Q = real(ifft(fGauss_Q,N)); % Correlated quadratures ...Gaussian vector
%% Normalizing the correlated gaussian vectors
tGauss_In = tGauss_I.*sqrt(1/(2*m*var(tGauss_I)));
tGauss_Qn = tGauss_Q.*sqrt(1/(2*m*var(tGauss_Q)));

end

