function x = Kalman_estimate(y_measured,v,ss_model_d,nt,nx,Q,R,P_prior)

Ad = ss_model_d.A;
Bd = ss_model_d.B;
Cd = ss_model_d.C;

% initialize x and P
x_prev = zeros(nx,1);
P = P_prior;
x  = zeros(nx,nt);

for i = 2:nt
    % PREDICTION STEP:
    %----------------------------------------------------------------------
    xPred    = Ad*x_prev + Bd*v(:,i-1);
    PPred    = Q + Ad*P*Ad';
    
    % CORRECTION STEP:
    %----------------------------------------------------------------------
    yPred    = y_measured(:,i-1)-Cd*xPred;
    S        = Cd*PPred*Cd' + R;
    K        = PPred*Cd'*inv(S);
    x_prev   = xPred + K*yPred;
    x(:,i)   = x_prev;
    P        = PPred - K*Cd*PPred;
end

end