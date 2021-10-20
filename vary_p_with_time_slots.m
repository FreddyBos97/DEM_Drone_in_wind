%% Caclulate the SSE of DEM for varying p to indicate the advantage
% of generalized coordinates figure 4 c
clear all
close all
clc

p_range = 0:6;

%% Main parameters
d_main      = 2;    % order of generalized coordinates for inputs
%N_AR        = 1;    % order of the AR system for SMIKF
N_SA        = 6;    % order of the AR system for SA
s_main      = 0.006;
Pz_main     = inv(8.1214e-09); % From determine noise for exp 25

% Settings for the input. For known inputs, sigma v is exp(-16)
sigma_v_main = diag([exp(-16) ones(1,3)*exp(-16)]);

% Settings for the time slots per experiment and the start and end times
N_slots = 5;
T_begin = 400;
T_end   = T_begin + 1200;
slots   = round(linspace(T_begin,T_end,N_slots+1));

exp_no_wind = [1,3,5,7];   % exp 21, exp 22 WM0 and exp 24 25 WM0
exp_wind    = [2,4,6,8]; % exp 21, exp 22 WM1, exp 24, exp 25 WM2

for j = 1:8
%% Load Data for phi and phi_dot, File numbers below
%                  21  22   24  25  26
% Wind mode    0   1   3     5   7  9
%              1   2   4     
%              2             6   8  10
% file 9 (exp26 WM0) is corrupted
    % shortest files are 1800
    % 

file_num = j;

for time_slot = 1:N_slots
    
Data = load_data(file_num,slots(time_slot),slots(time_slot+1)-1);

%% Convert the data to a model, containing the proper names and states
model = get_model_white_box(Data,0);

%% Find the proper noise charactaristics 
ms_num = 1; % number of multistarts for optimizing the s value
run_ms = 0; % choose 0 to skip the multistart 
model  = get_noise_charact(model,ms_num,run_ms);

% s is set to the sample time to provide a general method
model.s     = s_main;
model.d     = d_main; % Embedding of the inputs

model.sigma_v     = sigma_v_main;                % Very small for known input
model.prior_cause = model.v;                     % prior is the input for known causes
model.Pw          = model.Pw;                    % From the get_noise file
model.Pz          = Pz_main;             

%% Determine the SSE od DEM for varying p 
for i = 1:length(p_range)
    model.p  = p_range(i); % Embedding of the outputs
    brain = get_brain(model);
    
    %% DEM state estimation
    [out.x_DEM,model,brain] = DEM_Estimate(model,brain);
    
    
    %% Second order SA
    if p_range(i) == 0
        N_AR = 1;
    else
        N_AR = p_range(i);
    end 
    kalman.P_prior = eye(model.nx);
    kalman.Q = inv(model.Pw);
    kalman.R = inv(model.Pz);

    [SA2.AR_sigma,SA2.AR_par,SA2.AR_noise] = fit_AR(model.w,N_AR);
    [SA2.sys_aug_d,SA2.cov_w_aug] = augmented_kalman(N_AR,model.nx,...
        model.nv,model.ny,SA2.AR_par,model.sys_d,SA2.AR_sigma);
    SA2.P_prior_aug = eye(model.nx*(N_AR+1));
    out.x_SA2 = Kalman_estimate(model.y_meas,model.v,SA2.sys_aug_d,...
        model.nt,model.nx*(N_AR+1),SA2.cov_w_aug,kalman.R,SA2.P_prior_aug);
    
    %% Determine SSE for the hidden state phi dot
    SSE.trim = 10; % Trim of the inaccurate values at the edges
    SSE_DEM(time_slot,i)   = determine_sse(model.x_meas(2,:),out.x_DEM(2,:),SSE.trim);
    SSE_SA2(time_slot,i)   = determine_sse(model.x_meas(2,:),out.x_SA2(2,:),SSE.trim);
end
    SSE_main_DEM{j} =SSE_DEM;
    SSE_main_SA2{j} =SSE_SA2;
end
end

%% Results 
DEM_table = SSE_main_DEM;

for k = 1:8
%     for l = 1:length(p_range)
%         th = 2*median(DEM_table{k}(:,l));
%         DEM_table{k}(DEM_table{k}(:,l)>=th,l) = NaN
%     end
    SSE_DEM_mean(k,:) = mean(DEM_table{k},1,'omitnan');
    SSE_SA2_mean(k,:) = mean(SSE_main_SA2{k},1);
end

%% Plot the figures 
% plot figure with error bars

mean_no_wind = mean(SSE_DEM_mean(exp_no_wind,:));
mean_wind = mean(SSE_DEM_mean(exp_wind,:));
std_no_wind = std(SSE_DEM_mean(exp_no_wind,:));
std_wind = std(SSE_DEM_mean(exp_wind,:));

SA2_mean_wind = mean(SSE_SA2_mean(exp_wind,:));

x_plot = linspace(0,6.5,100);
g = fittype('a-b*exp(-c*x)');

expo_fit_no_wind = fit(p_range.',mean_no_wind.',g);
expo_fit_wind = fit(p_range.',mean_wind.',g);

fig3 = figure;
hold on 
xlim([0,6.5])
xlabel('Embedding order p')
ylabel('State estimation SSE')
ylim([0,25])

p1 = errorbar(p_range,mean_wind,std_wind,'Markersize',2);
p3 = plot(x_plot,expo_fit_wind(x_plot),'Color',[0, 0.4470, 0.7410],'LineWidth',2);
% p4 = plot(p_range,SA2_mean_wind,'x-');

set(p1, {'DisplayName'}, {'mean SSE with standard deviation'})
set(p3, {'DisplayName'}, {'Exponential fit'})
legend('Interpreter','latex')
ax = gca; % current axes
ax.FontSize = 15;

%% Save
saveas(fig3,'Figures/vary_p.eps','epsc')
saveas(fig3,'Figures/vary_p.jpg','jpg')
saveas(fig3,'Figures/vary_p.fig','fig')
