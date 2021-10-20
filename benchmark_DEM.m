%% Bechmark DEM against Kalman filter, state augmentation and first order 
% SMIKF, for one experiment. Use the DEM files from Ajith A. 
% Meera without altering the files. Only compare for the states phi and 
% phi dot Figure 4 a
clear all
close all
clc

%% Load Data for phi and phi_dot, File numbers below
%                  21  22   24  25  26
% Wind mode    0   1   3     5   7  9
%              1   2   4     
%              2             6   8  10
% file 9 (exp26 WM0) is corrupted
    % shortest files are 1800
    % 

file_num = 2;
T_begin = 400;
T_end   = T_begin+1200;
Data = load_data(file_num,T_begin,T_end);

%% Convert the data to a model, containing the proper names and states
model = get_model_white_box(Data,0);

%% Set up the properties for DEM, SA and SMIKF
model.p  = 6; % Embedding of the outputs
model.d  = 2; % Embedding of the inputs
model.N  = 1; % Order of the AR model for SA and SMIKF
model.N2 = 6; % Order of the higher order AR model for SA

%% Find the proper noise charactaristics 
ms_num = 1; % number of multistarts for optimizing the s value
run_ms = 0; % choose 0 to skip the multistart 
model  = get_noise_charact(model,ms_num,run_ms);
% s is set to the sample time to provide a general method
model.s = 0.006;


model.sigma_v     = eye(model.nv)*exp(-16); 
model.prior_cause = model.v;
model.Pw          = model.Pw;                    % From the get_noise file

model.Pz          = inv(8.1214e-09);             % From determine noise for exp 25
%model.Pz          = inv(9.83e-9);                % From Dennis Benders' thesis

brain = get_brain(model);
%% State estimation with DEM 
[out.x_DEM,model,brain] = DEM_Estimate(model,brain);
[Pi_xx,Pi_vv] = DEM_precision(brain);

%% State estimation Kalman
kalman.P_prior = eye(model.nx);
kalman.Q = inv(model.Pw);
kalman.R = inv(model.Pz);
out.x_KF = Kalman_estimate(model.y_meas,model.v,model.sys_d,model.nt,...
    model.nx,kalman.Q,kalman.R,kalman.P_prior);

%% Second order SA
[SA.AR_sigma,SA.AR_par,SA.AR_noise] = fit_AR(model.w,model.N);
[SA.sys_aug_d,SA.cov_w_aug] = augmented_kalman(model.N,model.nx,...
    model.nv,model.ny,SA.AR_par,model.sys_d,SA.AR_sigma);

[SA2.AR_sigma,SA2.AR_par,SA2.AR_noise] = fit_AR(model.w,model.N2);
[SA2.sys_aug_d,SA2.cov_w_aug] = augmented_kalman(model.N2,model.nx,...
    model.nv,model.ny,SA2.AR_par,model.sys_d,SA2.AR_sigma);
SA2.P_prior_aug = eye(model.nx*(model.N2+1));
out.x_SA2 = Kalman_estimate(model.y_meas,model.v,SA2.sys_aug_d,...
    model.nt,model.nx*(model.N2+1),SA2.cov_w_aug,kalman.R,SA2.P_prior_aug);

%% State estimation using SMIKF
SMIKF.P_prior_SMIKF{1} = eye(model.nx);
SMIKF.cov_w_SMIKF = diag(SA.AR_sigma.^2);
% for experiment 4 the AR parameters cuase the SMIKF algorithm to be unstable
% so for that experiment the AR parameters have been multiplied with 0.9
if file_num == 1
out.x_SMIKF = SMIKF1(model.y_meas,model.v,model.sys_d,model.nt,...
    model.nx,SMIKF.cov_w_SMIKF,kalman.R,SMIKF.P_prior_SMIKF,SA.AR_par*AR_mod);
else
    out.x_SMIKF = SMIKF1(model.y_meas,model.v,model.sys_d,model.nt,...
    model.nx,SMIKF.cov_w_SMIKF,kalman.R,SMIKF.P_prior_SMIKF,SA.AR_par);
end

%% Plot states
figure
hold on
plot(model.T,model.x_meas(1,:))
plot(model.T,out.x_DEM(1,:))
plot(model.T,out.x_KF(1,:))
plot(model.T,out.x_SA2(1,:))
plot(model.T,out.x_SMIKF(1,:))
title('$\phi$','Interpreter','latex')
legend('$\phi_{m}$','$\phi_{DEM}$','$\phi_{KF}$','$\phi_{SA}$',...
    '$\phi_{SMIKF}$','Interpreter','latex')

%% 
state_est_example = figure;
hold on
LW = 2;
plot(model.T,model.x_meas(2,:),'LineWidth',LW)
plot(model.T,out.x_DEM(2,:),'LineWidth',LW)
plot(model.T,out.x_KF(2,:),'--','LineWidth',1)
% plot(model.T,out.x_SA(2,:),'LineWidth',LW)
plot(model.T,out.x_SA2(2,:),'--','LineWidth',1)
plot(model.T,out.x_SMIKF(2,:),'--','LineWidth',1)
% title('$\dot \phi$','Interpreter','latex')
legend('$\dot \phi_{meas}$','$\dot \phi_{DEM}$','$\dot \phi_{KF}$',...
    '$\dot \phi_{SA}$','$\dot \phi_{SMIKF}$',...
    'Interpreter','latex')
ax = gca; % current axes
ax.FontSize = 15;
xlabel('$Time[s]$','Interpreter','latex')
ylabel('$\dot \phi[rad]$','Interpreter','latex')
xlim([2.5,5.5])
ylim([-0.8,1.1])
%% Plot estimate with precision 
sigma_phi_dot = sqrt(inv(Pi_xx(2,2)));

figure
hold on
LW = 2;
plot(model.T,model.x_meas(2,:),'LineWidth',LW)
plot(model.T,out.x_DEM(2,:),'LineWidth',LW)
plot(model.T,out.x_DEM(2,:)+sigma_phi_dot,'Color',[0.5 0.5 0.5])
plot(model.T,out.x_DEM(2,:)-sigma_phi_dot,'Color',[0.5 0.5 0.5])
title('$\dot \phi$','Interpreter','latex')
legend('$\dot \phi_{m}$','$\dot \phi_{DEM}$','Interpreter','latex')
ax = gca; % current axes
ax.FontSize = 30;

%% Determine SSE for the hidden state phi dot
SSE.trim = 5; % Trim of the inaccurate values at the edges
SSE.KF = determine_sse(model.x_meas(2,:),out.x_KF(2,:),SSE.trim);
SSE.DEM = determine_sse(model.x_meas(2,:),out.x_DEM(2,:),SSE.trim);
SSE.SA2 = determine_sse(model.x_meas(2,:),out.x_SA2(2,:),SSE.trim);
SSE.SMIKF = determine_sse(model.x_meas(2,:),out.x_SMIKF(2,:),SSE.trim);

display(["KF","DEM","SA2","SMIKF";SSE.KF,SSE.DEM,SSE.SA2,SSE.SMIKF])

%% save figure
saveas(state_est_example,'Figures/state_est_example.eps','epsc2')
saveas(state_est_example,'Figures/state_est_example.jpg','jpg')
saveas(state_est_example,'Figures/state_est_example.fig','fig')

