function [AR_sigma,AR_par,AR_noise] = fit_AR(w_meas,N)
% Fit an AR system to the nosie and output its characteristics

[nx,nt] = size(w_meas);
for k = 1:nx
    AR_sys = ar(w_meas(k,:),N);
    AR_sigma(k) = sqrt(AR_sys.NoiseVariance);
    AR_par(k,:) = -AR_sys.A(2:end);
    AR_w = AR_sigma(k)*randn(1,nt);
    AR_noise(k,:) = zeros(1,nt);
    AR_noise(k,1:N) = AR_w(1:N);
    for i = N:nt-1
        AR_noise(k,i+1) = flip(AR_par(k,:))*AR_noise(k,i+1-N:i).' + AR_w(i);
    end
end 