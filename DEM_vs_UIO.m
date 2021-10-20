%% Input estimation with the drone data 
% DEM vs UIO performance figure 7
clear all
close all
clc

%% Main parameters
p_main      = 6;    % order of generalized coordinates for states and outputs
d_main      = 2;    % order of generalized coordinates for inputs
s_main      = 0.006;
Pz_main     = inv(8.1214e-09); % From determine noise for exp 25

% Settings for the input. For known inputs, sigma v is exp(-16)
sigma_v_main = diag([exp(0) ones(1,3)*exp(-16)]);
v_est = 1; % determine which iput to estimate
prior_cause_05 = 0.5; % make the prior cause of the to-be-estimated input 0.5
% N_step = 3;
UIO_gamma = 0.5;
observable_system = 1;

% Settings for the time slots per experiment and the start and end times
T_begin = 400;
T_end   = T_begin + 1200;

for i = 1:8

file_num = i;
Data = load_data(file_num,T_begin,T_end);

%% Convert the data to a model, containing the proper names and states
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
if prior_cause_05   
    model.prior_cause(v_est,:) = ones(v_est,model.nt)*1;
end
model.Pz                    = eye(model.ny)*Pz_main;
model.Pw = eye(2)*exp(3);

brain = get_brain(model);
%% State estimation with DEM with measured input
[out.x_DEM,model,brain] = DEM_Estimate(model,brain);
[Pi_xx,Pi_vv,Pi_xv] = DEM_precision(brain);

%% UIO
[UIO_z, UIO_v] = UIO_estimator(model.sys_d,model.x_meas,...
    model.y_meas,model.nt,model.v,0,UIO_gamma,v_est,model.v);

%% SSE vlaues 
trim = 10;

SSE_DEM(i) = determine_sse(out.x_DEM(15,:),model.v(v_est,:),trim);
SSE_UIO(i) = determine_sse(UIO_v,model.v(v_est,:),trim);

end 

%% Results 
exp_wind = [2,4,6,8];

display([SSE_DEM;SSE_UIO])

%% Bar plot 
bar_UIO = figure;
hold on
bar_mat = [SSE_DEM(exp_wind);SSE_UIO(exp_wind)];
% X = categorical({'KF','DEM','SA','SA6','SMIKF'});
% X = reordercats(X,{'KF','DEM','SA','SA6','SMIKF'});
X = categorical({'exp 1','exp 2','exp 3','exp 4'});
X = reordercats(X,{'exp 1','exp 2','exp 3','exp 4'});
h = bar(X,bar_mat.');
ylabel('Input SSE','Interpreter','latex')
legend('DEM','UIO')
ylim([0,57])
ax = gca;
ax.FontSize = 15;

%% 
saveas(bar_UIO,'Figures/bar_UIO.eps','epsc')
saveas(bar_UIO,'Figures/bar_UIO.jpg','jpg')
saveas(bar_UIO,'Figures/bar_UIO.fig','fig')


