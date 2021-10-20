function [W0,V0y,V0v] = precision_matrix(s_brain,P_brain_w,...
                                         P_brain_z,sigma_v,p,d,nv,nx,ny) 
% Function to determine the generalized precision matrices from A.A. Meera
                                     
                                     
if p>d, pp = p; else pp = d; end

% Create the smoothness matrix
k          = 0:pp;
x          = sqrt(2)*s_brain;
r(1 + 2*k) = cumprod(1 - 2*k)./(x.^(2*k));
  
Cov     = [];
for i = 1:pp+1;
    Cov = [Cov; r([1:pp+1] + i - 1)];
    r = -r;
end

% Cov1 = equilibrate(Cov);

% [P,R,C] = equilibrate(Cov); 
% cond(Cov);
% cond(R*P*Cov*C)
% Cov1 = R*P*Cov*C;
Cov_inv = inv(Cov);

W0 = kron(Cov_inv(1:p+1,1:p+1),P_brain_w);
V0y = kron(Cov_inv(1:p+1,1:p+1),P_brain_z);
V0v = kron(Cov_inv(1:d+1,1:d+1),eye(nv)/sigma_v^2);

end