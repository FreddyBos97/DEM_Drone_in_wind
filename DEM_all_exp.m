%% Bechmark DEM against Kalman filter, state augmentation and first order 
% SMIKF, for all experiments. Use the DEM files from Ajith A. 
% Meera without altering the files. Only compare for the states phi and 
% phi dot Figure 4 b
clear all
close all
clc

%% Main parameters
p_main      = 6;    % order of generalized coordinates for states and outputs
d_main      = 2;    % order of generalized coordinates for inputs
N_AR        = 1;    % order of the AR system for SMIKF
N_SA        = 6;    % order of the AR system for SA
s_main      = 0.006;
Pz_main     = inv(8.1214e-09); % From determine noise for exp 25

% Settings for the input. For known inputs, sigma v is exp(-16)
sigma_v_main = diag([exp(-16) ones(1,3)*exp(-16)]);
v_est = 1; % determine which iput to estimate

% Settings for the time slots per experiment and the start and end times
N_slots = 5;
T_begin = 400;
T_end   = T_begin + 1200;
slots   = round(linspace(T_begin,T_end,N_slots+1));

no_wind = [1,3,5,7];   % exp 21, exp 22 WM0 and exp 24 25 WM0
wind    = [2,4,6,8]; % exp 21, exp 22 WM1, exp 24, exp 25 WM2
%% Loop through all the data and save the SSE and the SSE per second 
for experiment = 1:10
%% Load Data for phi and phi_dot, File numbers below
%                  21  22   24  25  26
% Wind mode    0   1   3     5   7  9
%              1   2   4     
%              2             6   8  10
% file 9 (exp26 WM0) is corrupted
    % shortest files are 1800 samples long

file_num = experiment;
for time_slot = 1:N_slots
    
Data = load_data(file_num,slots(time_slot),slots(time_slot+1)-1);
%% Convert the data to a model, containing the proper names and states
model = get_model_white_box(Data,0);

%% Set up the properties for DEM, SA and SMIKF
model.p  = p_main; % Embedding of the outputs
model.d  = d_main; % Embedding of the inputs
model.N  = N_AR;   % Order of the AR model for SA and SMIKF
model.N2 = N_SA;   % Order of the higher order AR model for SA

%% Find the proper noise charactaristics 
ms_num = 1; % number of multistarts for optimizing the s value
run_ms = 0; % choose 0 to skip the multistart 
model  = get_noise_charact(model,ms_num,run_ms);

% s is set to the sample time to provide a general method
model.s = s_main;
AR_mod = 1;

model.sigma_v     = sigma_v_main;                  % Very small for known input
model.prior_cause = model.v;
model.Pw          = model.Pw;                      % From the get_noise file
model.Pz          = eye(model.ny)*Pz_main;

brain = get_brain(model);
% brain.Pw = eye(2)*exp(-4);
%% State estimation with DEM 
[out.x_DEM,model,brain] = DEM_Estimate(model,brain);

%% State estimation Kalman
kalman.P_prior = eye(model.nx);
kalman.Q = inv(model.Pw);
kalman.R = inv(model.Pz);
out.x_KF = Kalman_estimate(model.y_meas,model.v,model.sys_d,model.nt,...
    model.nx,kalman.Q,kalman.R,kalman.P_prior);

%% State estimation State Augmentation (SA)
[SA.AR_sigma,SA.AR_par,SA.AR_noise] = fit_AR(model.w,model.N);
[SA.sys_aug_d,SA.cov_w_aug] = augmented_kalman(model.N,model.nx,...
    model.nv,model.ny,SA.AR_par,model.sys_d,SA.AR_sigma);

[SA2.AR_sigma,SA2.AR_par,SA2.AR_noise] = fit_AR(model.w,model.N2);

[SA2.sys_aug_d,SA2.cov_w_aug] = augmented_kalman(model.N2,model.nx,...
    model.nv,model.ny,SA2.AR_par*AR_mod,model.sys_d,SA2.AR_sigma);
SA2.P_prior_aug = eye(model.nx*(model.N2+1));

out.x_SA2 = Kalman_estimate(model.y_meas,model.v,SA2.sys_aug_d,...
    model.nt,model.nx*(model.N2+1),SA2.cov_w_aug,kalman.R,SA2.P_prior_aug);

%% State estimation using SMIKF
SMIKF.P_prior_SMIKF{1} = eye(model.nx);
SMIKF.cov_w_SMIKF = diag(SA.AR_sigma.^2);

out.x_SMIKF = SMIKF1(model.y_meas,model.v,model.sys_d,model.nt,...
    model.nx,SMIKF.cov_w_SMIKF,kalman.R,SMIKF.P_prior_SMIKF,SA.AR_par*AR_mod);

%% Determine SSE for the hidden state phi dot
SSE.trim    = 10; % Trim of the inaccurate values at the edges
SSE.KF      = determine_sse(model.x_meas(2,:),out.x_KF(2,:),SSE.trim);
SSE.DEM     = determine_sse(model.x_meas(2,:),out.x_DEM(2,:),SSE.trim);
SSE.SA2     = determine_sse(model.x_meas(2,:),out.x_SA2(2,:),SSE.trim);
SSE.SMIKF   = determine_sse(model.x_meas(2,:),out.x_SMIKF(2,:),SSE.trim);

% display(["KF","DEM","SA","SA2","SMIKF";SSE.KF,SSE.DEM,SSE.SA,SSE.SA2,SSE.SMIKF])

% store the data 
save_data{experiment,time_slot}.SSE = SSE;
save_data{experiment,time_slot}.model = model;
save_data{experiment,time_slot}.brain = brain;
save_data{experiment,time_slot}.SA2 = SA2;
save_data{experiment,time_slot}.SMIKF = SMIKF;
save_data{experiment,time_slot}.kalman = kalman;
save_data{experiment,time_slot}.out = out;
save_data{experiment,time_slot}.wind = Data.wind_vel;

estimated_v{experiment,time_slot} = out.x_DEM([(model.p+1)*model.nx+1:(model.p+1)*model.nx + model.nv],:);
end
end 

%keyboard
%% Examine data
%% Determine noise std

for i = 1:8
    for j = 1:N_slots
    std_phi(i,j)   = std(save_data{i,j}.model.x_meas(1,:));
    std_phid(i,j)  = std(save_data{i,j}.model.x_meas(2,:));
    std_w1(i,j)    = std(save_data{i,j}.model.w(1,:));
    std_w2(i,j)    = std(save_data{i,j}.model.w(2,:));
    
    m_phi(i,j)   = mean(save_data{i,j}.model.x_meas(1,:));
    m_phid(i,j)  = mean(save_data{i,j}.model.x_meas(2,:));
    m_w1(i,j)    = mean(save_data{i,j}.model.w(1,:));
    m_w2(i,j)    = mean(save_data{i,j}.model.w(2,:));
    end 
end

std_phi_nw  = mean(mean(std_phi(no_wind,:),2));
std_phid_nw = mean(mean(std_phid(no_wind,:),2));
std_w1_nw   = mean(mean(std_w1(no_wind,:),2));
std_w2_nw   = mean(mean(std_w2(no_wind,:),2));

std_phi_w  = mean(mean(std_phi(wind,:),2));
std_phid_w = mean(mean(std_phid(wind,:),2));
std_w1_w   = mean(mean(std_w1(wind,:),2));
std_w2_w   = mean(mean(std_w2(wind,:),2));

display([" ","phi","phi d","w1","w2";...
    "Without wind",std_phi_nw,std_phid_nw,std_w1_nw,std_w2_nw;...
    "With wind",std_phi_w,std_phid_w,std_w1_w,std_w2_w])

%% No wind data
clear SSE_no_wind_table
for i = 1:length(no_wind)
    for j = 1:N_slots
        k =  N_slots*(i-1) + j;
        SSE_no_wind_table(k,:) = cell2mat(struct2cell(save_data{no_wind(i),j}.SSE));
    end 
end
SSE_no_wind_table = SSE_no_wind_table(:,[2,3,4,5]);
SSE_no_wind_sum = mean(SSE_no_wind_table);
display(["KF","DEM","SA","SMIKF";SSE_no_wind_table;SSE_no_wind_sum])

%% Wind data
clear SSE_wind_table
for i = 1:length(no_wind)
    for j = 1:N_slots
        k =  N_slots*(i-1) + j;
        SSE_wind_table(k,:) = cell2mat(struct2cell(save_data{wind(i),j}.SSE));
    end 
end

SSE_wind_table = SSE_wind_table(:,[2,3,4,5]);
% Remove the largest outliers
out_thresh = 100;
SSE_wind_table(SSE_wind_table>=out_thresh) = NaN;

SSE_wind_sum = mean(SSE_wind_table,'omitnan');
display(["KF","DEM","SA","SMIKF";SSE_wind_table;SSE_wind_sum])

figure
hold on
plot(SSE_wind_table)
%% Make the bar plot

fig_bar = figure;
bar_mat = [SSE_wind_sum([1,4,3,2]);SSE_no_wind_sum([1,4,3,2])];
X = categorical({'With Wind','Without Wind'});
X = reordercats(X,{'With Wind','Without Wind'});
h = bar(X,bar_mat);
ax = gca; 
ax.FontSize = 15;
h(1).FaceColor = [0.9290, 0.6940, 0.1250];
h(2).FaceColor = [0.4660, 0.6740, 0.1880];
h(3).FaceColor = [0.4940, 0.1840, 0.5560];
h(4).FaceColor = [0.8500, 0.3250, 0.0980];

ylabel('State Estimation SSE')
set(h, {'DisplayName'}, {'KF','SMIKF','SA','DEM'}')
legend('Interpreter','latex') 
set(gca, 'YScale', 'log')
ylim([0.15 4])


%% save the data 
saveas(fig_bar,'Figures/benchmark_bar.eps','epsc')
saveas(fig_bar,'Figures/benchmark_bar.jpg','jpg')
saveas(fig_bar,'Figures/benchmark_bar.fig','fig')
