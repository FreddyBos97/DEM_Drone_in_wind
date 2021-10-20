%% Accuracy vs complexity, figures 8 and 9
clear all
close all
clc

%% Main parameters
p_main      = 6;    % order of generalized coordinates for states and outputs
d_main      = 2;    % order of generalized coordinates for inputs
s_main      = 0.006;
Pz_main     = inv(8.1214e-09); % From determine noise for exp 25
P_w_main    = eye(2)*exp(3);

v_est = 1; % determine which iput to estimate
prior_cause = 1;
UIO_gamma = 0.5;
observable_system = 0;

% Settings for the time slots per experiment and the start and end times
T_begin = 400;
T_end   = T_begin + 1200;

vary_w = linspace(-10 ,10 , 41);

exp_wind = [2,4,6,8];

for exp_num = 1:length(exp_wind)
    file_num = exp_wind(exp_num);
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
    
for i = 1:length(vary_w)
    P_v_1       = exp(vary_w(i));
    sigma_v_1 = 1/sqrt(P_v_1);
    sigma_v_2 = exp(-16);
    
    % Settings for the input. For known inputs, sigma v is exp(-16)
    sigma_v_main = diag([sigma_v_1 sigma_v_2 ones(1,2)*exp(-16)]);

    model.sigma_v       = sigma_v_main; 

    model.prior_cause           = model.v;
    model.prior_cause(v_est,:)  = ones(1,model.nt)*prior_cause;
    model.Pz                    = eye(model.ny)*Pz_main;
    model.Pw = P_w_main;

    brain = get_brain(model);
    %% State estimation with DEM 
    [out.x_DEM,model,brain] = DEM_Estimate(model,brain);
    [Pi_xx,Pi_vv,Pi_xv] = DEM_precision(brain);

    %% Determine SSE's 
    trim = 10;
    SSE_DEM_state(exp_num,i) = determine_sse(model.x_meas(2,:),out.x_DEM(2,:),trim);
    SSE_DEM_input(exp_num,i) = determine_sse(model.v(v_est,:),out.x_DEM(...
        (model.p+1)*model.nx+v_est,:),trim);

    v_DEM{exp_num}(i,:) = out.x_DEM((model.p+1)*model.nx+v_est,:);

    % Second input
        SSE_DEM_input2(exp_num,i) = determine_sse(model.v(2,:),out.x_DEM(...
        (model.p+1)*model.nx+2,:),trim);
        v_DEM2{exp_num}(i,:) = out.x_DEM((model.p+1)*model.nx+2,:);
        real_input2(exp_num,:) = model.v(2,:);
        
    x_DEM{exp_num}(i,:) = out.x_DEM(2,:);
    real_input(exp_num,:) = model.v(v_est,:);
    real_state(exp_num,:) = model.x_meas(2,:);

end 
    %% UIO
    [UIO_z, UIO_v(exp_num,:)] = UIO_estimator(model.sys_d,model.x_meas,...
        model.y_meas,model.nt,model.v,0,UIO_gamma,v_est,model.v);

    SSE_UIO_input(exp_num) = determine_sse(model.v(v_est,:),UIO_v(exp_num,:),trim);
end

%% Plot for in the paper
SSE_vs_Pv = figure;
hold on

p1 = plot(exp(vary_w),SSE_DEM_state,'LineWidth',2,'Color',...
    [0, 0.4470, 0.7410],'DisplayName','State estimation SSE');
p2 = plot(exp(vary_w),SSE_DEM_input,'LineWidth',2,'Color',...
    [0.8500, 0.3250, 0.0980],'DisplayName','DEM input estimation SSE');

h = [p1(1),p2(1)];
legend(h,'Interpreter','latex','Location','NorthWest')

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

%% Examine the input estimations

for k = 1:exp_num
figure
plot(model.T,real_input(k,:),'DisplayName','real input')
hold on
for j = 1:length(vary_w)
    plot(model.T,v_DEM{k}(j,:),'--','DisplayName',['sigma_v = ',num2str(vary_w(j))])
end
legend()
title(['Input prediction of exp ',num2str(k)])
end
plot(model.T,UIO_v,'DisplayName','UIO')

%% Show the effect of Pv on the states
Pv.Data = load_data(2,T_begin,T_end);
Pv.model = get_model_white_box(Pv.Data,0);
Pv.model.p  = p_main; % Embedding of the outputs
Pv.model.d  = d_main; % Embedding of the inputs
Pv.model  = get_noise_charact(Pv.model,1,0);

Pv.model.s             = s_main;
Pv.model.prior_cause           = Pv.model.v;
Pv.model.prior_cause(v_est,:) = ones(v_est,Pv.model.nt)*1;
Pv.model.sigma_v = sigma_v_main;
Pv.model.Pz                    = eye(Pv.model.ny)*Pz_main;
Pv.model.Pw = eye(2)*P_w_main;

% vary sigma v 
Pv.sigma_v_vary = [exp(-8),exp(-5),exp(-3),exp(0)]; 
Pv.P_vary = 1./Pv.sigma_v_vary.^2;

for j = 1:length(Pv.sigma_v_vary)
    Pv.model.sigma_v(v_est,v_est) = Pv.sigma_v_vary(j); 
    Pv.brain = get_brain(Pv.model);
    [Pv.x_DEM_vary_sigma_v{j},Pv.model,Pv.brain] = DEM_Estimate(Pv.model,Pv.brain);
    Pv.v_DEM(j,:) = Pv.x_DEM_vary_sigma_v{j}((Pv.model.p+1)*Pv.model.nx+v_est,:);
    Pv.real_input(j,:) = Pv.model.v(v_est,:);
end
%% Plot 
input_Pv = figure;
hold on
plot(Pv.model.T,Pv.model.v(v_est,:),'--','LineWidth',2)
plot(Pv.model.T,Pv.v_DEM,'LineWidth',2)
legend('Measured Input','$P_v = e^{16}$','$P_v = e^{10}$','$P_v = e^{6}$','$P_v = e^{0}$','Interpreter','latex')
ylim([-0.6,1.5])
xlabel('Time[s]')
ylabel('Input')
ax = gca;
ax.FontSize = 15;

%% Save figures
saveas(SSE_vs_Pv,'Figures/SSE_vs_Pv.eps','epsc2')
saveas(SSE_vs_Pv,'Figures/SSE_vs_Pv.jpg','jpg')
saveas(SSE_vs_Pv,'Figures/SSE_vs_Pv.fig','fig')

saveas(input_Pv,'Figures/input_Pv.eps','epsc2')
saveas(input_Pv,'Figures/input_Pv.jpg','jpg')
saveas(input_Pv,'Figures/input_Pv.fig','fig')
