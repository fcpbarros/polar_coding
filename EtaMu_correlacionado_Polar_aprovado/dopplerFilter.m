function [H] = dopplerFilter(fd, fs, N)

fdRatio = fd/fs; % Ratio between the sampling frequency and do ... Doppler shift [dimensionless]
km = floor(fdRatio*N);

H = zeros(1,N); % Filter's frequency responde

for k = 1:N
    if k == 1
        H(k) = 0;
    elseif k >= 2 && k <= km
        H(k) = sqrt(1./(2*sqrt(1-((k-1)/(N*fdRatio)).^2)));
    elseif k == km+1
        H(k) = sqrt((km/2)*((pi/2)-atan((km-1)/sqrt(2*km-1))));
    elseif k >= km+2 && k <= N-km
        H(k) = 0;
    elseif k == N-km+1
        H(k) = sqrt((km/2)*((pi/2)-atan((km-1)/sqrt(2*km-1))));
    else
        H(k) = sqrt(1./(2*sqrt(1-((N-(k-1))/(N*fdRatio)).^2)));
    end
end

end





