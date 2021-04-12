function etaMuChannel_Envelope = kappamu_corre_coefs(N,mu,kappa,Var,Lambda,fd,fs)

[etaMuChannel_I, etaMuChannel_Q] = kappaMuChannelGen(fd, fs, N, mu*2, ...
Var, Lambda, kappa);%, fJointComponents_Theory, xStep, ...
%markovSignEstimation);

% Composing the channel signal
etaMuChannel = etaMuChannel_I+1i*etaMuChannel_Q;
% Define the signal envelope
etaMuChannel_Envelope = abs(etaMuChannel');

end

