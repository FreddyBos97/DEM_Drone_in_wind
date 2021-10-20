%% Input estimation with the drone data, figure 6
clear all
close all
clc

%% Main parameters
p_main      = 6;    % order of generalized coordinates for states and outputs
d_main      = 2;    % order of generalized coordinates for inputs
s_main      = 0.006;
Pz_main     = inv(8.1214e-09); % From determine noise for exp 25
P_w_man     = exp(3);
manual_P_w  = 1;

% Settings for the input. For known inputs, sigma v is exp(-16)
sigma_v_main = diag([exp(0) ones(1,3)*exp(-16)]);
v_est = 1; % determine which iput to estimate
UIO_gamma = 0.5;
observable_system = 1;

% Settings for the time slots per experiment and the start and end times
T_begin = 400;
T_end   = T_begin + 1200;

file_num = 4;
Data = load_data(file_num,T_begin,T_end);

%% Convert the data to a model, containing the proper names and states
% use an input of 1 for fully observed states for UIO
model = get_model_white_box(Data,observable_system);

%% Set up the properties for DEM, SA and SMIKF
model.p  = p_main; % Embedding of the outputs
model.d  = d_main; % Embedding of the inputs


%% Find the proper noise charactaristics 
ms_num = 1; % number of multistarts for optimizing the s value
run_ms = 0; % choose 0 to skip the multistart 
model  = get_noise_charact(model,ms_num,run_ms);

model.s             = s_main;
model.sigma_v       = sigma_v_main; 

model.prior_cause           = model.v;
model.prior_cause(1,:)      = ones(v_est,model.nt)*0.5;
model.Pz                    = eye(model.ny)*Pz_main;

if manual_P_w
    model.Pw = eye(2)*P_w_man;
end

brain = get_brain(model);
%% State estimation with DEM 
[out.x_DEM,model,brain] = DEM_Estimate(model,brain);
[Pi_xx,Pi_vv,Pi_xv] = DEM_precision(brain);

%% UIO
[UIO_z, UIO_v] = UIO_estimator(model.sys_d,model.x_meas,...
    model.y_meas,model.nt,model.v,0,UIO_gamma,v_est,model.v);
%% plot input
fig_DEM_vs_UIO = figure;
hold on
plot(model.T,model.v(v_est,:),'LineWidth',2)
plot(model.T,out.x_DEM((model.p+1)*model.nx+v_est,:),'LineWidth',2)
plot(model.T,UIO_v,'--','LineWidth',2)
ax = gca; % current axes
ax.FontSize = 15;

ylabel('v')
xlabel('Time[s]')

legend('Measured v','DEM v','UIO v','Interpreter','latex')

%% Save figures
saveas(fig_DEM_vs_UIO,'Figures/input_estimation.eps','epsc2')
saveas(fig_DEM_vs_UIO,'Figures/input_estimation.jpg','jpg')
saveas(fig_DEM_vs_UIO,'Figures/input_estimation.fig','fig')

