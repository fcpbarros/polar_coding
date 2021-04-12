function coefs = alphamu_func(N,alpha,mu)

        %%Parametros
        omega=1; 
        rms=omega.^(1/(alpha));
        sigma2=(rms.^alpha)./(2.*mu);
        %% Geração dos coeficientes
        coefs = 0;
        for i = 1:2*mu
            coefs = coefs + (sqrt(sigma2).*randn(N,1)).^2;
        end
        coefs = coefs.^(1./alpha);

end

