function [etaMuChannel_I,etaMuChannel_Q]= etaMuChannelGen(fd, fs, ...
    N, m, Var, CorCoef, eta)%, fJointComponents_Theory, xStep, ...
    %markovSignEstimation)

    
%% m = 2*mu
H = dopplerFilter(fd, fs, N);

%% Generate the Eta-Mu signal
% Initializing variables
sumGauss_I = zeros(1,N);
sumGauss_Q = zeros(1,N);
tEtaMu_I = zeros(1,N);
tEtaMu_Q = zeros(1,N);

%% Generate the Covariance matriz
% Define the variance of the in-phase component
Var_I = 2*Var*eta/(1+eta);

% Define the variance of the quadrature component
Var_Q = 2*Var/(1+eta);

% Building the covariance matriz
matrizCov = [Var_I CorCoef*sqrt(Var_I)*sqrt(Var_Q) ; ...
    CorCoef*sqrt(Var_I)*sqrt(Var_Q) Var_Q];

% Performing the Cholesky decomposition for correlated variates
UpperChol = chol(matrizCov);


for cont = 1:m
    % Generate the Gaussian samples
    [tGauss_In, tGauss_Qn] = gaussGen(H, N, m);
    
    % Applying the correlation for the in-phase and in-quadrature ...components
    tGauss_Transposed = [tGauss_In' tGauss_Qn'];
    tGauss_Correlated = tGauss_Transposed*UpperChol;
    tGauss_In = tGauss_Correlated(:,1)';
    tGauss_Qn = tGauss_Correlated(:,2)';
    
    % The summation of the square gaussian signals compose the ...Eta-Mu signal
    tEtaMu_I = tEtaMu_I + (tGauss_In).^2;
    tEtaMu_Q = tEtaMu_Q + (tGauss_Qn).^2;
    
    % The summation of the signals defines the signal of the Eta-Mu
    % signal
    sumGauss_I = sumGauss_I + tGauss_In;
    sumGauss_Q = sumGauss_Q + tGauss_Qn;
end

%% Check if m is a non-integer
% If it's true, the in-phase and in-quadrature components will have ...a m+1
% cluster with a fraction of the samples with filled with zeros

if mod(m,1) ~= 0
    % Define the fraction of zeros samples
    Nmod = (1/mod(m,1));
    Ni = floor(N/Nmod);
    % Generate the Gaussian samples
    [tGauss_In, tGauss_Qn] = gaussGen(H, N, m);
    tGauss_In = [tGauss_In(1:Ni) zeros(1,N-Ni)];
    if m > 1
        tGauss_Qn = [tGauss_Qn(1:Ni) zeros(1,N-Ni)];
    else
        tGauss_Qn = [zeros(1,N-Ni) tGauss_Qn(N-Ni+1:N)];
    end
    
    % Applying the correlation for the in-phase and in-quadrature ...components
    tGauss_Transposed = [tGauss_In' tGauss_Qn'];
    tGauss_Correlated = tGauss_Transposed*UpperChol;
    tGauss_In = tGauss_Correlated(:,1)';
    tGauss_Qn = tGauss_Correlated(:,2)';
    % The summation of the square gaussian signals compose the ...Eta-Mu signal
    tEtaMu_I = tEtaMu_I + (tGauss_In).^2;
    tEtaMu_Q = tEtaMu_Q + (tGauss_Qn).^2;
    
    % The summation of the signals defines the signal of the Eta-Mu
    % signal
    sumGauss_I = sumGauss_I + tGauss_In;
    sumGauss_Q = sumGauss_Q + tGauss_Qn;
end

%% Define the signal estimation type
% if markovSignEstimation && (m) ~= 1
%     % Defining the signal of the Eta-Mu signal (Markov Chain)
%     [sign_I, sign_Q] = signalEstimationMarkov(N, ...
%         fJointComponents_Theory, xStep);
% else
%     % Defining the signal of the Eta-Mu signal
%     sign_I = sign(sumGauss_I);
%     sign_Q = sign(sumGauss_Q);
% end

sign_I = sign(sumGauss_I);
sign_Q = sign(sumGauss_Q);

tEtaMu_I = sign_I.*sqrt(tEtaMu_I);
tEtaMu_Q = sign_Q.*sqrt(tEtaMu_Q);

% Composing the channel signal
etaMuChannel_I = tEtaMu_I;
etaMuChannel_Q = tEtaMu_Q;

end

