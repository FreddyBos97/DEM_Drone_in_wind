%% Free Energy  Principle for State and Input Estimation of a 
%% QuadcopterFlying in Wind. A Matlab script to obtain the results and 
%% figures. 
%% Fred Bos, TU Delft, with use of the DEM code from Ajith Anil Meera and 
%% the data files from Dennis Benders.

%% List of results and figures in the paper
%%      Table:      1
%%      Figures:    3(a,b,c), 4(a,b,c), 5, 6, 7, 8, 9

%% The data folder contains 10 files with data: experiments 21, 22, 24, 25 
%% and 26 with both no wind(NW) and with wind(WW) conditions.
close all
clear all
clc

%% Part 1: Noise characterisation(NC)
% figures 3(a,b,c)

% select the start and end time in samples
NC.start_time  = 400; 
NC.end_time    = 800;
% select the correct file numbers
NC.Data_files_NW = [1,3,5,7];   % exp 21, 22, 24, 25 NW
NC.Data_files_WW = [2,4,6,8];   % exp 21, 22, 24, 25 WW

for i = 1:4
    % Select the correct file
    NC.data_file_NW{i} = NC.Data_files_NW(i);
    NC.data_file_WW{i} = NC.Data_files_WW(i);
    % Load the data
    NC.Data_NW{i} = load_data(NC.data_file_NW{i},NC.start_time,NC.end_time);
    NC.Data_WW{i} = load_data(NC.data_file_WW{i},NC.start_time,NC.end_time); 
    % Convert the data to the model format that is used in all scripts
    NC.Data_NW{i} = get_model_white_box(NC.Data_NW{i},0);
    NC.Data_WW{i} = get_model_white_box(NC.Data_WW{i},0);
    % Calculate the noise vectors
    NC.Data_NW{i} = get_noise_charact(NC.Data_NW{i},0,0);
    NC.Data_WW{i} = get_noise_charact(NC.Data_WW{i},0,0);

    % Check the standard deviations of the nosie
    NC.std_NW_w_phi(i)    = std(NC.Data_NW{i}.w(1,:));
    NC.std_NW_w_phid(i)   = std(NC.Data_NW{i}.w(2,:));
    NC.std_WW_w_phi(i)    = std(NC.Data_WW{i}.w(1,:));
    NC.std_WW_w_phid(i)   = std(NC.Data_WW{i}.w(2,:));

    NC.std_NW_phi(i)    = std(NC.Data_NW{i}.x_meas(1,:));
    NC.std_NW_phid(i)   = std(NC.Data_NW{i}.x_meas(2,:));
    NC.std_WW_phi(i)    = std(NC.Data_WW{i}.x_meas(1,:));
    NC.std_WW_phid(i)   = std(NC.Data_WW{i}.x_meas(2,:));
end

%% Plot the figures and obtain the table
NC.std_table = [mean(NC.std_NW_phi),mean(NC.std_NW_phid),...
                mean(NC.std_NW_w_phi),mean(NC.std_NW_w_phid);...
                mean(NC.std_WW_phi),mean(NC.std_WW_phid),...
                mean(NC.std_WW_w_phi),mean(NC.std_WW_w_phid)];
disp([[" ","phi","phi_d","w_phi","w_phi_d"];...
    [["Without Wind";"With wind"],NC.std_table]]);

%% Plot the Autocorrelation of the noise for with wind experiments 
NC.Plot_colors = [0, 0.4470, 0.7410   % Set the first four colours 
    0.8500, 0.3250, 0.0980
    0.9290, 0.6940, 0.1250
    0.4940, 0.1840, 0.5560];

NC.auto_lags = 50; % Determine the amount of lags for the autocorrelation
NC.font_size = 15; % Set the fontsize

auto_corr_WW = figure;

for l = 1:4 % Loop through the data to contstruct the autocorrelation plots
    subplot(2,1,1)
    hold on

    NC.a_c_phi(l,:) = autocorr(NC.Data_WW{l}.w(1,:),NC.auto_lags);

    plot(linspace(0,NC.auto_lags,NC.auto_lags+1),NC.a_c_phi(l,:),'.-',...
        'MarkerSize',NC.font_size)
    plot([linspace(0,NC.auto_lags,NC.auto_lags+1);...
        linspace(0,NC.auto_lags,NC.auto_lags+1)],...
        [zeros(1,NC.auto_lags+1);NC.a_c_phi(l,:)],...
        'Color',NC.Plot_colors(l,:))

    ylabel('$w_{\phi}$', 'Interpreter','latex') % Set the labels for the 
    title(' ')                                  % top plot
    NC.ax = gca; 
    NC.ax.FontSize = NC.font_size;

    subplot(2,1,2)
    hold on

    NC.a_c_phid(l,:) = autocorr(NC.Data_WW{l}.w(2,:),NC.auto_lags);
    NC.p(l) = plot(linspace(0,NC.auto_lags,NC.auto_lags+1),...
        NC.a_c_phid(l,:),'.-','MarkerSize',NC.font_size);
    plot([linspace(0,NC.auto_lags,NC.auto_lags+1);...
        linspace(0,NC.auto_lags,NC.auto_lags+1)],...
        [zeros(1,NC.auto_lags+1);NC.a_c_phid(l,:)],...
        'Color',NC.Plot_colors(l,:))
end

% Set the legend of the plot and set the label for the bottom plot
legend([NC.p(1),NC.p(2),NC.p(3),NC.p(4)],...
    {'Wind Exp 1','Wind Exp 2','Wind Exp 3','Wind Exp 4'})
NC.main_legend = legend([NC.p(1),NC.p(2),NC.p(3),NC.p(4)],...
    {'Wind Exp 1','Wind Exp 2','Wind Exp 3','Wind Exp 4'});
NC.newPosition = [0.65 0.25 0.25 0.25];
NC.newUnits = 'normalized';
set(NC.main_legend,'Position',NC.newPosition,'Units',NC.newUnits);

ylabel('$w_{\dot{\phi}}$', 'Interpreter','latex')
title(' ')
NC.ax = gca; 
NC.ax.FontSize = NC.font_size;

%% Construct the histogram plots 
NC.hist_bins = 15; 

hist_NW = figure; % Histogram with Gaussian fit for no wind condition 
hold on

subplot(1,2,1)
hold on
histfit(NC.Data_NW{2}.w(1,:),NC.hist_bins);
xlabel('$w_{\phi}$', 'Interpreter','latex')
NC.ax = gca; 
NC.ax.FontSize = NC.font_size;

subplot(1,2,2)
hold on
histfit(NC.Data_NW{2}.w(2,:),NC.hist_bins);
xlabel('$w_{\dot{\phi}}$', 'Interpreter','latex')
NC.ax = gca; 
NC.ax.FontSize = NC.font_size;

hist_WW = figure; % Histogram with Gaussian fit for wind condition
hold on

subplot(1,2,1)
hold on
histfit(NC.Data_WW{2}.w(1,:),NC.hist_bins);
xlabel('$w_{\phi}$', 'Interpreter','latex')
xlim([-2e-3,2e-3])
NC.ax = gca; 
NC.ax.FontSize = NC.font_size;

subplot(1,2,2)
hold on
histfit(NC.Data_WW{2}.w(2,:),NC.hist_bins);
xlabel('$w_{\dot{\phi}}$', 'Interpreter','latex')
xlim([-0.25,0.25])
NC.ax = gca; 
NC.ax.FontSize = NC.font_size;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: State estimation 
% figures 4(a,b,c), 5
% Main parameters
clear l i

SE.p_main      = 6;     % order of generalized coordinates for states and outputs
SE.d_main      = 2;     % order of generalized coordinates for inputs
SE.N_AR        = 1;     % order of the AR system for SMIKF
SE.N_SA        = 6;     % order of the AR system for SA
SE.s_main      = 0.006; % Smoothness value for DEM 
SE.Pz_main     = inv(8.1214e-09); % From determine noise for exp 25

% Settings for the input. For known inputs, sigma v is exp(-16)
SE.sigma_v_main = diag(ones(1,4)*exp(-16));

% Settings for the time slots per experiment and the start and end times
% Divide the experiments in smaller timeslots 
SE.N_slots = 5;
SE.T_begin = 400;
SE.T_end   = SE.T_begin + 1200;
SE.slots   = round(linspace(SE.T_begin,SE.T_end,SE.N_slots+1));
SSE.trim    = 10; % Trim of the inaccurate values at the edges
vary_p.trim = 10; % Trim of the inaccurate values at the edges
    
SE.no_wind = [1,3,5,7]; % exp 21, exp 22 WM0 and exp 24 25 WM0
SE.wind    = [2,4,6,8]; % exp 21, exp 22 WM1, exp 24, exp 25 WM2

% settings for varying embedding order p 
vary_p.p_range = 0:6;

% Loop through all the data and save the SSE and the SSE per second 
for experiment = 1:10
%% Load Data for phi and phi_dot, File numbers below
%                  21  22   24  25  26
% Wind mode    0   1   3     5   7  9
%              1   2   4     
%              2             6   8  10
% file 9 (exp26 WM0) is corrupted
    % shortest files are 1800 samples long

SE.file_num = experiment;
for time_slot = 1:SE.N_slots
    
SE.Data_slots = load_data(SE.file_num,SE.slots(time_slot),SE.slots(time_slot+1)-1);
% Convert the data to a model, containing the proper names and states
model = get_model_white_box(SE.Data_slots,0);

% Set up the properties for DEM, SA and SMIKF
model.p  = SE.p_main; % Embedding of the outputs
model.d  = SE.d_main; % Embedding of the inputs
model.N  = SE.N_AR;   % Order of the AR model for SA and SMIKF
model.N2 = SE.N_SA;   % Order of the higher order AR model for SA
model.s  = SE.s_main; % Smoothness value for the Gaussian noise 

% Find the proper noise charactaristics 
SE.ms_num = 1; % number of multistarts for optimizing the s value
SE.run_ms = 0; % choose 0 to skip the multistart 
model  = get_noise_charact(model,SE.ms_num,SE.run_ms);

model.sigma_v     = SE.sigma_v_main;          % Very small for known input
model.prior_cause = model.v;
model.Pw          = model.Pw;                 % From the get_noise file
model.Pz          = eye(model.ny)*SE.Pz_main;

brain = get_brain(model); % Setup the brain component for DEM 

% State estimation with DEM 
[out.x_DEM,model,brain] = DEM_Estimate(model,brain);

% State estimation with varying p
vary_p.model = model;
for i = 1:length(vary_p.p_range)
    vary_p.model.p  = vary_p.p_range(i); % Embedding of the outputs
    vary_p.brain = get_brain(vary_p.model);
    
    % DEM state estimation
    [vary_p.x_DEM,vary_p.model,vary_p.brain] = DEM_Estimate(vary_p.model,vary_p.brain);
     
    % Determine SSE for the hidden state phi dot
    vary_p.SSE_DEM(time_slot,i)   = determine_sse(vary_p.model.x_meas(2,:)...
        ,vary_p.x_DEM(2,:),vary_p.trim);
end

% State estimation Kalman
kalman.P_prior = eye(model.nx);
kalman.Q = inv(model.Pw);
kalman.R = inv(model.Pz);
out.x_KF = Kalman_estimate(model.y_meas,model.v,model.sys_d,model.nt,...
    model.nx,kalman.Q,kalman.R,kalman.P_prior);

% %% State estimation State Augmentation (SA)
% [SA.AR_sigma,SA.AR_par,SA.AR_noise] = fit_AR(model.w,model.N);
% [SA.sys_aug_d,SA.cov_w_aug] = augmented_kalman(model.N,model.nx,...
%     model.nv,model.ny,SA.AR_par,model.sys_d,SA.AR_sigma);
% 
% SA.P_prior_aug = eye(model.nx*(model.N+1));
% out.x_SA = Kalman_estimate(model.y_meas,model.v,SA.sys_aug_d,model.nt,...
%     model.nx*(model.N+1),SA.cov_w_aug,kalman.R,SA.P_prior_aug);

% Second order SA
[SA2.AR_sigma,SA2.AR_par,SA2.AR_noise] = fit_AR(model.w,model.N2);
[SA2.sys_aug_d,SA2.cov_w_aug] = augmented_kalman(model.N2,model.nx,...
    model.nv,model.ny,SA2.AR_par,model.sys_d,SA2.AR_sigma);
SA2.P_prior_aug = eye(model.nx*(model.N2+1));
out.x_SA2 = Kalman_estimate(model.y_meas,model.v,SA2.sys_aug_d,...
    model.nt,model.nx*(model.N2+1),SA2.cov_w_aug,kalman.R,SA2.P_prior_aug);

% State estimation using SMIKF
[SMIKF.AR_sigma,SMIKF.AR_par,SMIKF.AR_noise] = fit_AR(model.w,model.N);
SMIKF.P_prior_SMIKF{1} = eye(model.nx);
SMIKF.cov_w_SMIKF = diag(SMIKF.AR_sigma.^2);
out.x_SMIKF = SMIKF1(model.y_meas,model.v,model.sys_d,model.nt,...
    model.nx,SMIKF.cov_w_SMIKF,kalman.R,SMIKF.P_prior_SMIKF,SMIKF.AR_par);

% Determine SSE for the hidden state phi dot
SSE.KF      = determine_sse(model.x_meas(2,:),out.x_KF(2,:),SSE.trim);
SSE.DEM     = determine_sse(model.x_meas(2,:),out.x_DEM(2,:),SSE.trim);
%SSE.SA      = determine_sse(model.x_meas(2,:),out.x_SA(2,:),SSE.trim);
SSE.SA2     = determine_sse(model.x_meas(2,:),out.x_SA2(2,:),SSE.trim);
SSE.SMIKF   = determine_sse(model.x_meas(2,:),out.x_SMIKF(2,:),SSE.trim);
% display(["KF","DEM","SA","SA2","SMIKF";SSE.KF,SSE.DEM,SSE.SA,SSE.SA2,SSE.SMIKF])

% store the data 
save_data{experiment,time_slot}.SSE = SSE;
save_data{experiment,time_slot}.model = model;
save_data{experiment,time_slot}.brain = brain;
%save_data{experiment,time_slot}.SA = SA;
save_data{experiment,time_slot}.SA2 = SA2;
save_data{experiment,time_slot}.SMIKF = SMIKF;
save_data{experiment,time_slot}.kalman = kalman;
save_data{experiment,time_slot}.out = out;
%save_data{experiment,time_slot}.wind = Data_slots.wind_vel;
%estimated_v{experiment,time_slot} = out.x_DEM([(model.p+1)*model.nx+1:(model.p+1)*model.nx + model.nv],:);
end
vary_p.all_SSE{experiment} = vary_p.SSE_DEM;
end 
% clear all the data that is just for one experiment for one time slot 
clear time_slot experiment SA2 SMIKF kalman brain model SSE out

%% Load data for the visualisation of state estimation 
Vis.Data = load_data(2,400,1200);
Vis.model = get_model_white_box(Vis.Data,0);

% Set up the properties for DEM, SA and SMIKF
Vis.model.p  = SE.p_main; % Embedding of the outputs
Vis.model.d  = SE.d_main; % Embedding of the inputs
Vis.model.N  = SE.N_AR;   % Order of the AR model for SA and SMIKF
Vis.model.N2 = SE.N_SA;   % Order of the higher order AR model for SA

% Find the proper noise charactaristics 
Vis.model   = get_noise_charact(Vis.model,1,0);
Vis.model.s = SE.s_main;

Vis.model.sigma_v     = SE.sigma_v_main; 
Vis.model.prior_cause = Vis.model.v;
Vis.model.Pz          = inv(8.1214e-09); % From determine noise for exp 25

Vis.brain = get_brain(Vis.model);
% State estimation with DEM 
[Vis.x_DEM,Vis.model,Vis.brain] = DEM_Estimate(Vis.model,Vis.brain);

% State estimation Kalman
Vis.kalman.P_prior = eye(Vis.model.nx);
Vis.kalman.Q = inv(Vis.model.Pw);
Vis.kalman.R = inv(Vis.model.Pz);
Vis.x_KF = Kalman_estimate(Vis.model.y_meas,Vis.model.v,Vis.model.sys_d,...
    Vis.model.nt,Vis.model.nx,Vis.kalman.Q,Vis.kalman.R,Vis.kalman.P_prior);

% Second order SA
[Vis.SA2.AR_sigma,Vis.SA2.AR_par,Vis.SA2.AR_noise] =...
    fit_AR(Vis.model.w,Vis.model.N2);
[Vis.SA2.sys_aug_d,Vis.SA2.cov_w_aug] = augmented_kalman(Vis.model.N2,...
    Vis.model.nx,Vis.model.nv,Vis.model.ny,Vis.SA2.AR_par,Vis.model.sys_d,...
    Vis.SA2.AR_sigma);
Vis.SA2.P_prior_aug = eye(Vis.model.nx*(Vis.model.N2+1));
Vis.x_SA2 = Kalman_estimate(Vis.model.y_meas,Vis.model.v,Vis.SA2.sys_aug_d,...
    Vis.model.nt,Vis.model.nx*(Vis.model.N2+1),Vis.SA2.cov_w_aug,...
    Vis.kalman.R,Vis.SA2.P_prior_aug);

% State estimation using SMIKF
[Vis.SMIKF.AR_sigma,Vis.SMIKF.AR_par,Vis.SMIKF.AR_noise]...
    = fit_AR(Vis.model.w,Vis.model.N);
Vis.SMIKF.P_prior_SMIKF{1} = eye(Vis.model.nx);
Vis.SMIKF.cov_w_SMIKF = diag(Vis.SMIKF.AR_sigma.^2);
Vis.x_SMIKF = SMIKF1(Vis.model.y_meas,Vis.model.v,Vis.model.sys_d,...
    Vis.model.nt,Vis.model.nx,Vis.SMIKF.cov_w_SMIKF,Vis.kalman.R,...
    Vis.SMIKF.P_prior_SMIKF,Vis.SMIKF.AR_par);

%% Examine the SSE data 
SE.SSE_no_wind_table = zeros(SE.N_slots*length(SE.no_wind),5);
SE.SSE_wind_table    = zeros(SE.N_slots*length(SE.no_wind),5);

for i = 1:length(SE.no_wind)
    for j = 1:SE.N_slots
        k =  SE.N_slots*(i-1) + j;
        
        SE.SSE_no_wind_table(k,:) = cell2mat(struct2cell(...
            save_data{SE.no_wind(i),j}.SSE)); % no wind data
        SE.SSE_wind_table(k,:) = cell2mat(struct2cell(...
            save_data{SE.wind(i),j}.SSE)); % wind data
    end 
end
SE.SSE_no_wind_table = SE.SSE_no_wind_table(:,[2,3,4,5]);
SE.SSE_no_wind_sum = mean(SE.SSE_no_wind_table);
display(["KF","DEM","SA","SMIKF";SE.SSE_no_wind_table;SE.SSE_no_wind_sum])

SE.SSE_wind_table = SE.SSE_wind_table(:,[2,3,4,5]);
SE.out_thresh = 100; % Remove the largest outliers
SE.SSE_wind_table(SE.SSE_wind_table>=SE.out_thresh) = NaN;
SE.SSE_wind_sum = mean(SE.SSE_wind_table,'omitnan');
display(["KF","DEM","SA","SMIKF";SE.SSE_wind_table;SE.SSE_wind_sum])

%% Demonstrate DEM climbes the VFE curve
VFE.model = save_data{2,1}.model;
VFE.brain = save_data{2,1}.brain;
VFE.x_DEM = save_data{2,1}.out.x_DEM;
VFE.X_varied = linspace(-2,2,101);

VFE.x_tilde = VFE.x_DEM([1:(VFE.model.p+1)*VFE.model.nx],:);
VFE.y_tilde = VFE.model.Y_embed([1:(VFE.model.p+1)*VFE.model.ny],:);
VFE.v_tilde =  -VFE.model.Y_embed([1+(VFE.model.p+1)*VFE.model.ny:end],:);
VFE.eta_tilde = -VFE.model.Y_embed([1+(VFE.model.p+1)*VFE.model.ny:end],:);

VFE.eps_y = VFE.y_tilde - VFE.brain.Ct*VFE.x_tilde;
VFE.eps_v = VFE.v_tilde - VFE.eta_tilde;
VFE.eps_x = VFE.brain.Da*VFE.x_tilde - VFE.brain.At*VFE.x_tilde - VFE.brain.Bt*VFE.v_tilde;

VFE.eps_tilde = [VFE.eps_y;VFE.eps_v;VFE.eps_x];

VFE.Pi_tilde = [VFE.brain.V0y,zeros(size(VFE.brain.V0y,1),size(VFE.brain.V0v,1)),...
    zeros(size(VFE.brain.V0y,1),size(VFE.brain.W0,1));zeros(size(VFE.brain.V0v,1),...
    size(VFE.brain.V0y,1)),VFE.brain.V0v,zeros(size(VFE.brain.V0v,1),...
    size(VFE.brain.W0,1));zeros(size(VFE.brain.W0,1),size(VFE.brain.V0y,1)),...
    zeros(size(VFE.brain.W0,1),size(VFE.brain.V0v,1)),VFE.brain.W0];

for j = 1:VFE.model.nt
    VFE.V_real(j) = -0.5*VFE.eps_tilde(:,j).'*VFE.Pi_tilde*VFE.eps_tilde(:,j);
end

for k = 1:length(VFE.model.T)-1
for j = 1:length(VFE.X_varied)
VFE.x_meas_var = VFE.model.x_meas;
VFE.x_meas_var(2,k) = VFE.X_varied(j);

VFE.x_meas_embed = zeros((VFE.model.p+1)*VFE.model.nx,VFE.model.nt);

for i = 1:length(VFE.model.T)
    if i>VFE.model.p+1 && i<length(VFE.model.T)-VFE.model.p-1
        VFE.x_meas_embed(:,i) = embed_Y(VFE.x_meas_var,VFE.model.p+1,VFE.model.T(i),VFE.model.sam_time);
    else
        VFE.x_meas_embed([1,2],i) = VFE.x_meas_var(:,i);
        VFE.x_meas_embed([3:end],i) = 0;
    end
end

VFE.x_tilde1 = VFE.x_meas_embed(:,k+1);
VFE.y_tilde1 = VFE.y_tilde(:,k+1);
VFE.v_tilde1 = VFE.eta_tilde(:,k+1);
VFE.eta_tilde1 = VFE.eta_tilde(:,k+1);

VFE.eps_y = VFE.y_tilde1 - VFE.brain.Ct*VFE.x_tilde1;
VFE.eps_v = VFE.v_tilde1 - VFE.eta_tilde1;
VFE.eps_x = VFE.brain.Da*VFE.x_tilde1 - VFE.brain.At*VFE.x_tilde1 - VFE.brain.Bt*VFE.v_tilde1;
VFE.eps_tilde = [VFE.eps_y;VFE.eps_v;VFE.eps_x];

VFE.V(k+1,j) = -0.5*VFE.eps_tilde.'*VFE.Pi_tilde*VFE.eps_tilde;
end
end

%% Plot the visualisation of state estimation for DEM and the benchmarks
state_est_example = figure;
hold on
plot(Vis.model.T,Vis.model.x_meas(2,:),'LineWidth',2)
plot(Vis.model.T,Vis.x_DEM(2,:),'LineWidth',2)
plot(Vis.model.T,Vis.x_KF(2,:),'--','LineWidth',1)
plot(Vis.model.T,Vis.x_SA2(2,:),'--','LineWidth',1)
plot(Vis.model.T,Vis.x_SMIKF(2,:),'--','LineWidth',1)
% title('$\dot \phi$','Interpreter','latex')
legend('$\dot \phi_{meas}$','$\dot \phi_{DEM}$','$\dot \phi_{KF}$',...
    '$\dot \phi_{SA}$','$\dot \phi_{SMIKF}$',...
    'Interpreter','latex')
Vis.ax = gca; % current axes
Vis.ax.FontSize = 15;
xlabel('$Time[s]$','Interpreter','latex')
ylabel('$\dot \phi[rad]$','Interpreter','latex')
xlim([2.5,5.5])
ylim([-0.8,1.1])

% Make the bar plot
fig_bar = figure;
SE.bar_mat = [SE.SSE_wind_sum([1,4,3,2]);SE.SSE_no_wind_sum([1,4,3,2])];
SE.X = categorical({'With Wind','Without Wind'});
SE.X = reordercats(SE.X,{'With Wind','Without Wind'});
SE.h = bar(SE.X,SE.bar_mat);
SE.ax = gca; 
SE.ax.FontSize = 15;
SE.h(1).FaceColor = [0.9290, 0.6940, 0.1250];
SE.h(2).FaceColor = [0.4660, 0.6740, 0.1880];
SE.h(3).FaceColor = [0.4940, 0.1840, 0.5560];
SE.h(4).FaceColor = [0.8500, 0.3250, 0.0980];

ylabel('State Estimation SSE')
set(SE.h, {'DisplayName'}, {'KF','SMIKF','SA','DEM'}')
legend('Interpreter','latex') 
set(gca, 'YScale', 'log')
ylim([0.15 4])

%% Plot the embedding order vs the SSE
vary_p.DEM_table = vary_p.all_SSE;
for l = 1:8
    vary_p.SSE_DEM_mean(l,:) = mean(vary_p.DEM_table{l},1,'omitnan');
end

% plot figure with error bars
vary_p.mean_no_wind = mean(vary_p.SSE_DEM_mean(SE.no_wind,:));
vary_p.mean_wind = mean(vary_p.SSE_DEM_mean(SE.wind,:));
vary_p.std_no_wind = std(vary_p.SSE_DEM_mean(SE.no_wind,:));
vary_p.std_wind = std(vary_p.SSE_DEM_mean(SE.wind,:));

vary_p.x_plot = linspace(0,6.5,100);
vary_p.g = fittype('a-b*exp(-c*x)');

vary_p.expo_fit_no_wind = fit(vary_p.p_range.',vary_p.mean_no_wind.',vary_p.g);
vary_p.expo_fit_wind = fit(vary_p.p_range.',vary_p.mean_wind.',vary_p.g);

Vary_p_figure = figure;
hold on 
xlim([0,6.5])
xlabel('Embedding order p')
ylabel('State estimation SSE')
ylim([0,25])

vary_p.p1 = errorbar(vary_p.p_range,vary_p.mean_wind,vary_p.std_wind,'Markersize',2);
vary_p.p3 = plot(vary_p.x_plot,vary_p.expo_fit_wind(vary_p.x_plot),...
    'Color',[0, 0.4470, 0.7410],'LineWidth',2);

set(vary_p.p1, {'DisplayName'}, {'mean SSE with standard deviation'})
set(vary_p.p3, {'DisplayName'}, {'Exponential fit'})
legend('Interpreter','latex')
vary_p.ax = gca; % current axes
vary_p.ax.FontSize = 15;

%% VFE surface plot 
VFE.trim = 20;
VFE.cut_off_range = VFE.trim:VFE.model.nt-VFE.trim-40;

Free_energy_curve = figure;
hold on
VFE.V_surf = surf(VFE.X_varied,VFE.model.T(VFE.cut_off_range),VFE.V(VFE.cut_off_range,:));
VFE.V_surf.EdgeAlpha = 0.1;
VFE.p_DEM = plot3(VFE.x_DEM(2,VFE.cut_off_range),VFE.model.T(VFE.cut_off_range),...
    VFE.V_real(VFE.cut_off_range),'LineWidth',2);
legend([VFE.V_surf,VFE.p_DEM],{'Free Energy Curves','DEM State Estimate'}...
    ,'interpreter','latex','Location','NorthWest')

view(45,50)
caxis(1.0e+08*[-2 -0.0001])
VFE.ax = gca;
VFE.ax.FontSize = 15;

zlabel('VFE','Interpreter','latex')
xlabel('Varied $\dot \phi$[rad]','interpreter','latex')
ylabel('Time[s]','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Input estimation 
% Figures 6, 7, 8, 9
clear j k i l

% Main parameters
IE.p_main      = SE.p_main;    
IE.d_main      = SE.d_main; 
IE.s_main      = SE.s_main;
IE.Pz_main     = SE.Pz_main; 
IE.Pw_main     = eye(2)*exp(3);

% Settings for the input. For known inputs, sigma v is exp(-16) for known
% input 1 sigma is exp(0)
IE.sigma_v_main = diag([exp(0) ones(1,3)*exp(-16)]);
IE.P_exponent = linspace(-10 ,10 , 41);
IE.P_exponent = [IE.P_exponent,16]; % append one value for the plot
IE.v_est = 1; % determine which iput to estimate
IE.prior_cause_bench = 0.5; % make the prior cause of the to-be-estimated input 0.5
IE.prior_cause_acc = 1;     % for accuracy vs complexity the prior is 1 for a better visualisation
IE.UIO_gamma = 0.5;

% Settings for the time slots per experiment and the start and end times
IE.T_begin = 400;
IE.T_end   = IE.T_begin + 1200;
IE.trim    = 10;
IE.no_wind = [1,3,5,7]; % exp 21, exp 22 WM0 and exp 24 25 WM0
IE.wind    = [2,4,6,8]; % exp 21, exp 22 WM1, exp 24, exp 25 WM2

for i = 1:8

IE.file_num = i;
IE.Data = load_data(IE.file_num,IE.T_begin,IE.T_end);

% Convert the data to a model, containing the proper names and states
% In order to benchmark with UIO the system has to be fully observable
IE_obs.model = get_model_white_box(IE.Data,1);

% Set up the properties for DEM
IE_obs.model.p  = IE.p_main; % Embedding of the outputs
IE_obs.model.d  = IE.d_main; % Embedding of the inputs

% Find the proper noise charactaristics 
IE_obs.model  = get_noise_charact(IE_obs.model,1,0);

IE_obs.model.s             = IE.s_main;
IE_obs.model.sigma_v       = IE.sigma_v_main; 

IE_obs.model.prior_cause           = IE_obs.model.v;
IE_obs.model.prior_cause(1,:)      = ones(IE.v_est,IE_obs.model.nt)...
    *IE.prior_cause_bench;

IE_obs.model.Pz                    = eye(IE_obs.model.ny)*IE.Pz_main;
IE_obs.model.Pw = IE.Pw_main;

IE_obs.brain = get_brain(IE_obs.model);
% State estimation with DEM with measured input
[IE_obs.x_DEM,IE_obs.model,IE_obs.brain] = DEM_Estimate(IE_obs.model,...
    IE_obs.brain);

% UIO
[IE_obs.UIO_z, IE_obs.UIO_v] = UIO_estimator(IE_obs.model.sys_d,...
    IE_obs.model.x_meas,IE_obs.model.y_meas,IE_obs.model.nt,...
    IE_obs.model.v,0,IE.UIO_gamma,IE.v_est,IE_obs.model.v);

% SSE vlaues
IE_obs.SSE_DEM(i) = determine_sse(IE_obs.x_DEM(IE_obs.model.nx*...
    (IE_obs.model.p+1)+IE.v_est,:),IE_obs.model.v(IE.v_est,:),IE.trim);
IE_obs.SSE_UIO(i) = determine_sse(IE_obs.UIO_v,IE_obs.model.v(IE.v_est,:),IE.trim);

% save the data for plotting
IE.save_data{i}.model = IE_obs.model;
IE.save_data{i}.brain = IE_obs.brain;
IE.save_data{i}.x_DEM = IE_obs.x_DEM;
IE.save_data{i}.UIO_v = IE_obs.UIO_v;

% In order to show the accuracy vs complexity tradeoff the system cannot be
% fully observable
IE_uno.model = get_model_white_box(IE.Data,0); 
IE_uno.model.p  = SE.p_main; % Embedding of the outputs
IE_uno.model.d  = SE.d_main; % Embedding of the inputs
IE_uno.model.s  = SE.s_main;

for j = 1:length(IE.P_exponent)
IE_uno.P_v_1       = exp(IE.P_exponent(j));
IE_uno.sigma_v_1 = 1/sqrt(IE_uno.P_v_1);
IE_uno.model.sigma_v = diag([IE_uno.sigma_v_1 ones(1,3)*exp(-16)]); 

IE_uno.model.prior_cause          = IE_uno.model.v;
IE_uno.model.prior_cause(IE.v_est,:)...
    = ones(1,IE_uno.model.nt)*IE.prior_cause_acc;
IE_uno.model.Pz                   = eye(IE_uno.model.ny)*IE.Pz_main;
IE_uno.model.Pw                   = IE.Pw_main;

IE_uno.brain = get_brain(IE_uno.model);
% State estimation with DEM
[IE_uno.x_DEM,IE_uno.model,IE_uno.brain] = DEM_Estimate(IE_uno.model,IE_uno.brain);

% Determine SSE's 
IE_uno.SSE_DEM_state(i,j) = determine_sse(IE_uno.model.x_meas(2,:),...
    IE_uno.x_DEM(2,:),IE.trim);
IE_uno.SSE_DEM_input(i,j) = determine_sse(IE_uno.model.v(IE.v_est,:),...
    IE_uno.x_DEM((IE_uno.model.p+1)*IE_uno.model.nx+IE.v_est,:),IE.trim);

IE_uno.v_DEM{i}(j,:) = IE_uno.x_DEM((IE_uno.model.p+1)...
    *IE_uno.model.nx+IE.v_est,:);
end 
IE_uno.real_input(i,:) = IE_uno.model.v(IE.v_est,:);
end 
%% Show the effect of Pv on the states
Pv.Data = load_data(4,IE.T_begin,IE.T_end);
Pv.model = get_model_white_box(Pv.Data,0);
Pv.model.p  = IE.p_main; % Embedding of the outputs
Pv.model.d  = IE.d_main; % Embedding of the inputs
Pv.model  = get_noise_charact(Pv.model,1,0);

Pv.model.s             = IE.s_main;
Pv.model.prior_cause           = Pv.model.v;
Pv.model.prior_cause(IE.v_est,:) = ones(IE.v_est,Pv.model.nt)*IE.prior_cause_acc;
Pv.model.sigma_v = IE.sigma_v_main;
Pv.model.Pz                    = eye(Pv.model.ny)*IE.Pz_main;
Pv.model.Pw = eye(2)*IE.Pw_main;

% vary sigma v 
Pv.sigma_v_vary = [exp(-8),exp(-5),exp(-3),exp(0)]; 
Pv.P_vary = 1./Pv.sigma_v_vary.^2;

for j = 1:length(Pv.sigma_v_vary)
    Pv.model.sigma_v(IE.v_est,IE.v_est) = Pv.sigma_v_vary(j); 
    Pv.brain = get_brain(Pv.model);
    [Pv.x_DEM_vary_sigma_v{j},Pv.model,Pv.brain] = DEM_Estimate(Pv.model,Pv.brain);
    Pv.v_DEM(j,:) = Pv.x_DEM_vary_sigma_v{j}((Pv.model.p+1)*Pv.model.nx+IE.v_est,:);
    Pv.real_input(j,:) = Pv.model.v(IE.v_est,:);
end

%% Plot the UIO state and the DEM state
fig_DEM_vs_UIO = figure;
IE.plot_num = 2;
hold on
plot(IE.save_data{IE.plot_num}.model.T,IE.save_data{IE.plot_num}.model.v...
    (IE.v_est,:),'LineWidth',2)
plot(IE.save_data{IE.plot_num}.model.T,IE.save_data{IE.plot_num}.x_DEM...
    ((IE.save_data{IE.plot_num}.model.p+1)*IE.save_data{IE.plot_num}.model.nx...
    +IE.v_est,:),'LineWidth',2)
plot(IE.save_data{IE.plot_num}.model.T,IE.save_data{IE.plot_num}.UIO_v,...
    '--','LineWidth',2)
IE_obs.ax = gca; % current axes
IE_obs.ax.FontSize = 15;
ylabel('v')
xlabel('Time[s]')
legend('Measured v','DEM v','UIO v','Interpreter','latex')

%% Bar plot of DEM vs UIO 
bar_UIO = figure;
hold on
IE_obs.bar_mat = [IE_obs.SSE_DEM(IE.wind);IE_obs.SSE_UIO(IE.wind)];
IE_obs.X = categorical({'exp 1','exp 2','exp 3','exp 4'});
IE_obs.X = reordercats(IE_obs.X,{'exp 1','exp 2','exp 3','exp 4'});
IE_obs.h = bar(IE_obs.X,IE_obs.bar_mat.');
ylabel('Input SSE','Interpreter','latex')
legend('DEM','UIO')
IE_obs.ax = gca;
IE_obs.ax.FontSize = 15;
ylim([0 58])

%% Show the effect of Pv
input_Pv = figure;
hold on
plot(Pv.model.T,Pv.model.v(IE.v_est,:),'--','LineWidth',2)
plot(Pv.model.T,Pv.v_DEM,'LineWidth',2)
legend('Measured Input','$P_v = e^{16}$','$P_v = e^{10}$','$P_v = e^{6}$','$P_v = e^{0}$','Interpreter','latex')
ylim([-0.6,1.5])
xlabel('Time[s]')
ylabel('Input')
IE.ax = gca;
IE.ax.FontSize = 15;

%% Accuracy vs complexity plot 
SSE_vs_Pv = figure;
hold on

IE.p1 = plot(exp(IE.P_exponent(1:41)),IE_uno.SSE_DEM_state(IE.wind,[1:41]),...
    'LineWidth',2,'Color',[0, 0.4470, 0.7410],'DisplayName',...
    'State estimation SSE');
IE.p2 = plot(exp(IE.P_exponent(1:41)),IE_uno.SSE_DEM_input(IE.wind,[1:41]),...
    'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980],'DisplayName',...
    'DEM input estimation SSE');

IE.h = [IE.p1(1),IE.p2(1)];
legend(IE.h,'Interpreter','latex','Location','NorthWest')

ax1 = gca; % current axes

set(ax1, 'XScale', 'log')
set(ax1, 'YScale', 'log')
ylabel(ax1,'SSE','Interpreter','latex')
xlabel(ax1,'$P_{v1}$','Interpreter','latex')

ax1.FontSize = 15;
ax1_pos = ax1.Position; 

ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right', ...
    'Color','none') ;

ax2.XLim = ax1.XLim ;
ax2.YLim = ax1.YLim ;
set(ax2, 'XScale', 'log')
set(ax2, 'YScale', 'log')
ax2.FontSize = 15;

yticks(ax2,[5,230])
yticklabels(ax2,{'High Accuracy','Low Accuracy'})
ytickangle(ax2,90)
set(ax2,'TickLength',[0 0])

xticks(ax2,[0.001,1000])
xticklabels(ax2,{'Low Complexity','High Complexity'})
