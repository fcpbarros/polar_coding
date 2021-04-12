function [alphaMuChannel_I,alphaMuChannel_Q]= alphaMuChannelGen(fd, fs, ...
    N, m, Var, CorCoef, alpha)%, fJointComponents_Theory, xStep, ...
    %markovSignEstimation)

    
%% m = 2*mu
H = dopplerFilter(fd, fs, N);
omega=1; 
rms=omega.^(1/(alpha));
sigma2=(Var.^alpha)./(2.*m);

%% Generate the Eta-Mu signal
% Initializing variables
sumGauss_I = zeros(1,N);
sumGauss_Q = zeros(1,N);
talphaMu = zeros(1,N);
tEtaMu_Q = zeros(1,N);

%% Generate the Covariance matriz
% Define the variance of the in-phase component
Var_I = 2*Var*sigma2; %%testar Var_I = Var*sigma2

% 
% Define the variance of the quadrature component
Var_Q = 2*Var*sigma2; %%testar Var_Q = Var*sigma2

% Building the covariance matriz
matrizCov = [Var_I CorCoef*Var_I; ...
    CorCoef*Var_I Var_Q];

%matrizCov = nearestSPD(matrizCov);
% Performing the Cholesky decomposition for correlated variates
UpperChol = chol(matrizCov);


for cont = 1:m
    % Generate the Gaussian samples
    [tGauss_In, tGauss_Qn] = gaussGen(H, N, m,alpha);
    
    % Applying the correlation for the in-phase and in-quadrature ...components
    tGauss_Transposed = [tGauss_In' tGauss_Qn'];
    tGauss_Correlated = tGauss_Transposed*UpperChol;
    tGauss_In = tGauss_Correlated(:,1)';
    tGauss_Qn = tGauss_Correlated(:,2)';
    
    % The summation of the square gaussian signals compose the ...Eta-Mu signal
    talphaMu = talphaMu + (tGauss_In).^2;
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
    talphaMu = talphaMu + (tGauss_In).^2;
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

talphaMu = sign_I.*sqrt(talphaMu);
tEtaMu_Q = sign_Q.*sqrt(tEtaMu_Q);

% Composing the channel signal
alphaMuChannel_I = talphaMu;
alphaMuChannel_Q = tEtaMu_Q;

end

