function alphaMuChannel_Envelope = alphamu_corre_coefs(N,mu,alpha,Var,Lambda,fd,fs)

[alphaMuChannel_I,alphaMuChannel_Q]= alphaMuChannelGen(fd, fs, N, mu*2, ...
Var, Lambda, alpha);%, fJointComponents_Theory, xStep, ...
%markovSignEstimation);

% Composing the channel signal
alphaMuChannel = alphaMuChannel_I.^2 + alphaMuChannel_Q.^2;
% Define the signal envelope
alphaMuChannel_Envelope = ((alphaMuChannel')).^(1/alpha);

end

