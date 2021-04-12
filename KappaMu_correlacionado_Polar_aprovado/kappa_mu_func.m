function coefs = kappa_mu_func(N,kappa,mu)
    
%%Parametros
rms = 1;
sigma2=(rms.^2)./(2.*(kappa+1).*mu); 
P=sqrt(kappa.*sigma2);
%% Geração dos coeficientes
coefs = 0;
for i = 1:2*mu
    coefs = coefs + (sqrt(sigma2).*randn(N,1)+P).^2;
%     ytrans=sqrt(((sqrt(sigma2).*randn(1,n)+P).^2+(sqrt(sigma2)...
% .*randn(1,n)+P).^2)+(sqrt(sigma2).*randn(1,n)+P).^2);
end
coefs = sqrt(coefs);
end

