function etaMuChannel_Envelope = etamu_corre_coefs(N,mu,eta,Var,Lambda,fd,fs)

[etaMuChannel_I, etaMuChannel_Q] = etaMuChannelGen(fd, fs, N, mu*2, ...
Var, Lambda, eta);%, fJointComponents_Theory, xStep, ...
%markovSignEstimation);

% Composing the channel signal
etaMuChannel = etaMuChannel_I+1i*etaMuChannel_Q;
% Define the signal envelope
etaMuChannel_Envelope = abs(etaMuChannel');

end

