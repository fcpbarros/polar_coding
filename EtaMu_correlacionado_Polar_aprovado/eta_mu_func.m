function coefs = eta_mu_func(N,eta,mu)

%% parametros
rms = 1; 
sigmay2 =(rms.^2)./(2.*mu.*(1+eta));
sigmax2=sigmay2.*eta;
%% Geração dos coeficientes
coefs = 0;
for i = 1:2*mu
    coefs = coefs + (sqrt(sigmax2).*randn(N,1)).^2+(sqrt(sigmay2).*randn(N,1)).^2;
end
coefs = sqrt(coefs);
end

