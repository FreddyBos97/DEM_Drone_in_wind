% Funtion that handles the embedding of Y, from A.A. meera
% EDIT Generate generalized y - central approach: 2 front and 2 back
function dY = embed_Y(Y,n,t,dt)

    [~, N]  = size(Y);
    T = zeros(n);
    
    s      = round((t)/dt);
    k      = (1:n)  + fix(s - (n + 1)/2);
    x      = s - min(k) + 1;
    y_pad = k<x-1 | k> N-(x-1);% Extra step: truncate innacurate derivative at edges
    i      = k < 1;
    k      = k.*~i + i;
    i      = k > N;
    k      = k.*~i + i*N;


    % Inverse embedding operator (T): cf, Taylor expansion Y(t) <- T*y[:]
    %----------------------------------------------------------------------
    for i = 1:n
        for j = 1:n
            T(i,j) = ((i - x)*dt)^(j - 1)/prod(1:(j - 1));
        end
    end

    % embedding operator: y[:] <- E*Y(t)
    %----------------------------------------------------------------------
    E     = inv(T);
   
    dY      = Y(:,k)*E';
    dY(:,end-sum(y_pad)+1:end) = zeros(size(Y,1),sum(y_pad)); % truncate inaccurate derivatives
    dY = reshape(dY,[],1);
end