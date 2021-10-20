function [z_f, d_f] = UIO_estimator(sys,x_real,y_real,k,d,if_dataset,UIO_gamma,unobs_input,real_cause)

A = sys.A; C = sys.C; D = sys.D;
if size(sys.B,2) == 1
    B1 = zeros(size(A,1),1); B2 = sys.B;
else 
    B1 = sys.B(:,setdiff(1:end,unobs_input)); B2 = sys.B(:,unobs_input);
end

nx = size(A,1); ny = size(C,1); nu = size(B2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 1, as per the paper, solve the LMI directly and tune gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if if_dataset == 1
    P = sdpvar(nx,nx);
    Q = sdpvar(nu,nu);
    M = sdpvar(nu,nu);
    G = sdpvar(nx,ny);
    % gamma = 100;
    
    AA = [-P         zeros(nx,1)             (P*A-G*C)'     -(M*B2'*C'*C)'; ...
        zeros(1,nx)   -Q                      B2'*(P*A-G*C)' Q-B2'*(M*B2'*C'*C)'; ...
        (P*A-G*C)     (B2'*(P*A-G*C)')'       -P             zeros(nx,1); ...
        -(M*B2'*C'*C) Q'-(B2'*(M*B2'*C'*C)')' zeros(1,nx)    -Q];
    
    F = [P>=0, Q>=0, M>=0, AA<=0];
    optimize(F);
    P = value(P); Q = value(Q);
    M = value(M);
    G = value(G);
    L = inv(P)*G;
    gamma = M*inv(Q);
    L = value(L);
    gamma = value(gamma)*1200000; % tune for best gamma; motor = gamma*1200000

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Method 2, as per paper, stable observer design through eigen values
    %          tune parameter UIO_gamma for best performance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    L = sdpvar(nx,ny);
    gamma = sdpvar(nu,nu);
    AAA = [A-L*C (A-L*C)*B2; -gamma*B2'*C'*C eye(nu,nu)...
        - gamma*B2'*C'*C*B2];
    % optimize([sum(abs(eig(A-L*C)))<=1.75, ...
    %     sum(abs(eig(eye(ny)-C*B2*gamma*B2'*C')))<=.3, gamma>=.01, L>=0])
    optimize([sum(abs(eig(AAA)))<=size(AAA,1), gamma>=UIO_gamma])
    % optimize([max(abs(eig(AAA)))<=1, gamma>=UIO_gamma])
    
    L = value(L)
    gamma = value(gamma)   
    % value(abs(eig(AAA)))
end

if size(sys.B,2) == 1
    u = zeros(1,k);
else
    u = real_cause(setdiff(1:end,unobs_input),:);
end

z_real = x_real;

d_f = zeros(nu,k);
x_f = zeros(nx,k);
y_f = zeros(ny,k);
z_f = zeros(nx,k);

d_f0 = .01;
d_f(:,1) = d_f0 + gamma*B2'*C'*(y_real(:,1) - C*B2*d_f0); % x_f(1) = 0
x_f(:,2) = A*B2*d_f0 + B1*u(:,k) + L*(y_real(:,1) - C*B2*d_f0);

for i = 2:k-1
%     y_f(:,i) = C*x_f(:,i) + D*u(:,i) + C*B2*d_f(:,i-1);
    y_f(:,i) = C*x_f(:,i) + C*B2*d_f(:,i-1);
    d_f(:,i) = d_f(:,i-1) + gamma*B2'*C'*(y_real(:,i) - y_f(:,i));
    z_f(:,i) = x_f(:,i) + B2*d_f(:,i-1);
    x_f(:,i+1) = A*x_f(:,i) + B1*u(:,i) + A*B2*d_f(:,i-1) + ...
                                     L*(y_real(:,i) - y_f(:,i));
end
% y_f(:,k) = C*x_f(:,k) + D*u(:,k) + C*B2*d_f(:,k-1);
y_f(:,k) = C*x_f(:,k) + C*B2*d_f(:,k-1);
d_f(:,k) = d_f(:,k-1) + gamma*B2'*C'*(y_real(:,k) - y_f(:,k));
z_f(:,k) = x_f(:,k) + B2*d_f(:,k-1);

d_f(:,1:k-1) = d_f(:,2:k);
z_f(:,1:k-1) = z_f(:,2:k);
end