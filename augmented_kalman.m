function [ss_model_aug_d,cov_w_aug] = augmented_kalman(N,nx,nv,ny,AR_par,sys_d,AR_sigma)
% This functions extends the system equations with the AR parameters of the
% nosie for the state augmentation approach 

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;

nnx = N*nx;
AR_mat_top = zeros(nx,nnx);
for j = 1:N
    AR_mat_top(1:nx,1+(j-1)*nx:j*nx) = diag(AR_par(:,j));
end

AR_mat = [AR_mat_top;eye(nnx-nx,nnx)];

A_aug = [Ad,eye(nx,nnx);zeros(nnx,nx),AR_mat];
B_aug = [Bd;zeros(nnx,nv)];
C_aug = [Cd,zeros(ny,nnx)];

cov_w_aug = zeros(nx+nnx,nx+nnx);
cov_w_aug(nx+1:nx+nx,nx+1:nx+nx) = diag(AR_sigma.^2);

ss_model_aug_d = ss(A_aug,B_aug,C_aug,[]);

