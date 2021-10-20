function x_smikf = SMIKF1(y,v,ss_model_d,nt,nx,cov_w,cov_z,P_prior_SMIKF,AR_par)
% Function that performs the second moment information kalman filter
% state estimation

Ad = ss_model_d.A;
Bd = ss_model_d.B;
Cd = ss_model_d.C;

%% Initiate covariances and kalman gain
P{1} = P_prior_SMIKF{1};
Q = cov_w;
R = cov_z;
K = eye(nx,size(Cd,1));

% construct state vectors
x_smikf = zeros(nx,nt);
x_prior = zeros(nx,nt);

% Initiate covariance between noise and state
P_ww{1} = cov_w;
P_wx = P_prior_SMIKF{1};
P_xw = P_prior_SMIKF{1};

x_prev = zeros(nx,1);
phi = diag(AR_par);
for k = 2:nt
    % Prediction Step:
    xPred    = Ad*x_prev + Bd*v(:,k-1);
    
    P_ww{k} = phi^2 *  P_ww{k-1} + Q;
    P_w1w = phi*P_ww{k-1};
    
    P_xw = (eye(nx) - K*Cd)*P_w1w; 
    P_wx = P_xw.'; 
   
    PPred    = Ad*P{k-1}*Ad' + Ad*P_xw + P_wx*Ad.' + P_ww{k};
    
    % Correction Step:
    yPred    = y(:,k-1)-Cd*xPred;
    S        = Cd*PPred*Cd' + R;
    K        = PPred*Cd'*inv(S);
    x_prev   = xPred + K*yPred;
    
    x_smikf(:,k)   = x_prev;
    P{k}        = PPred - K*Cd*PPred;
end