clear; close; clc;
%% Desvanecimento Correlacionado Eta Mu
fd = 5;
%fd = 30;
fs = 60; % Sampling frequency [Hz]
N = 1024; % Number of samples of the Eta-Mu signal [dimensionless]
mu = 2;
%Var = 1; % Variance of the Eta-Mu signal / 0 < Var < Infinity
Lambda = 0.8;
kappa = 2;
Var=1;
[etaMuChannel_I, etaMuChannel_Q] = kappaMuChannelGen(fd, fs, N, mu*2, ...
Var, Lambda, kappa);%, fJointComponents_Theory, xStep, ...
%markovSignEstimation);
% Composing the channel signal
etaMuChannel = etaMuChannel_I+1i*etaMuChannel_Q;
% Define the signal envelope
etaMuChannel_Envelope = abs(etaMuChannel);
semilogy(etaMuChannel_Envelope);
legenda = sprintf('KappaMu %.2f/%.2f corr. fd %.2f ',kappa,mu,fd);
legend(legenda)
%save('2048etamu_fd60.mat','etaMuChannel_Envelope')
% Define the signal phase
% etaMuChannel_Phase = atan2(etaMuChannel_Q,etaMuChannel_I);