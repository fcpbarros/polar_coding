clear; close; clc;
%% Desvanecimento Correlacionado Eta Mu
fd = 5;
%fd = 30;
fs = 60; % Sampling frequency [Hz]
N = 512; % Number of samples of the Eta-Mu signal [dimensionless]
mu = 1;
Var = 1; % Variance of the Eta-Mu signal / 0 < Var < Infinity
Lambda = 0.8;
alpha = 2; 
[alphaMuChannel_I,alphaMuChannel_Q]= alphaMuChannelGen(fd, fs, N, mu*2, ...
Var, Lambda, alpha);%, fJointComponents_Theory, xStep, ...
%markovSignEstimation);

% Composing the channel signal
alphaMuChannel = alphaMuChannel_I.^2 + alphaMuChannel_Q.^2;
% Define the signal envelope
alphaMuChannel_Envelope = ((alphaMuChannel')).^(1/alpha);
semilogy(alphaMuChannel_Envelope);
legenda = sprintf('AlphaMu %.2f/%.2f corr. fd %.2f ',alpha,mu,fd);
legend(legenda)
%save('2048etamu_fd60.mat','etaMuChannel_Envelope')
% Define the signal phase
% etaMuChannel_Phase = atan2(etaMuChannel_Q,etaMuChannel_I);