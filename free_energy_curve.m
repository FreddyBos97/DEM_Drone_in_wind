%% DEM climbs the free energy curve, figure 5
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
T_end   = T_begin+200;
Data = load_data(file_num,T_begin,T_end);

%% Convert the data to a model, containing the proper names and states
model = get_model_white_box(Data,0);

%% Set up the properties for DEM, SA and SMIKF
model.p  = 6; % Embedding of the outputs
model.d  = 2; % Embedding of the inputs
trim     = 10;% Remove inaccurate edges

%% Find the proper noise charactaristics 
ms_num = 1; % number of multistarts for optimizing the s value
run_ms = 0; % choose 0 to skip the multistart 
model  = get_noise_charact(model,ms_num,run_ms);
% s is set to the sample time to provide a general method
%model.s = model.sam_time;
model.s = 0.006;

model.sigma_v     = eye(model.nv)*exp(-16); 
model.prior_cause = model.v;
model.Pw          = model.Pw;                    % From the get_noise file

model.Pz          = inv(8.1214e-09);             % From determine noise for exp 25
%model.Pz          = inv(9.83e-9);                % From Dennis Benders' thesis

brain = get_brain(model);
%% State estimation with DEM 
[x_DEM,model,brain] = DEM_Estimate(model,brain);
[Pi_xx,Pi_vv] = DEM_precision(brain);

figure
hold on
plot(model.T,model.x_meas(2,:))
plot(model.T,x_DEM(2,:))

%% Calculate the variational free energy of the estimated states

x_tilde = x_DEM([1:(model.p+1)*model.nx],:);
y_tilde = model.Y_embed([1:(model.p+1)*model.ny],:);
v_tilde =  -model.Y_embed([1+(model.p+1)*model.ny:end],:);
eta_tilde = -model.Y_embed([1+(model.p+1)*model.ny:end],:);

eps_y = y_tilde - brain.Ct*x_tilde;
eps_v = v_tilde - eta_tilde;
eps_x = brain.Da*x_tilde - brain.At*x_tilde - brain.Bt*v_tilde;

eps_tilde = [eps_y;eps_v;eps_x];
Pi_tilde = [brain.V0y,zeros(size(brain.V0y,1),size(brain.V0v,1)),...
    zeros(size(brain.V0y,1),size(brain.W0,1));zeros(size(brain.V0v,1),...
    size(brain.V0y,1)),brain.V0v,zeros(size(brain.V0v,1),...
    size(brain.W0,1));zeros(size(brain.W0,1),size(brain.V0y,1)),...
    zeros(size(brain.W0,1),size(brain.V0v,1)),brain.W0];
for j = 1:model.nt
    V_real(j) = -0.5*eps_tilde(:,j).'*Pi_tilde*eps_tilde(:,j);
end

V_real([1:trim,end-trim:end]) = 0;
figure
plot(V_real)

%% Vary one state and calculate the free energy 
clear V
X_varied = linspace(-2,2,101);


for k = 1:length(model.T)-1
for j = 1:length(X_varied)
x_meas_var = model.x_meas;
x_meas_var(2,k) = X_varied(j);

x_meas_embed = zeros((model.p+1)*model.nx,model.nt);

for i = 1:length(model.T)
    if i>model.p+1 && i<length(model.T)-model.p-1
        x_meas_embed(:,i) = embed_Y(x_meas_var,model.p+1,model.T(i),model.sam_time);
    else
        x_meas_embed([1,2],i) = x_meas_var(:,i);
        x_meas_embed([3:end],i) = 0;
    end
end

x_tilde1 = x_meas_embed(:,k+1);
y_tilde1 = y_tilde(:,k+1);
v_tilde1 = eta_tilde(:,k+1);
eta_tilde1 = eta_tilde(:,k+1);

eps_y = y_tilde1 - brain.Ct*x_tilde1;
eps_v = v_tilde1 - eta_tilde1;
eps_x = brain.Da*x_tilde1 - brain.At*x_tilde1 - brain.Bt*v_tilde1;
eps_tilde = [eps_y;eps_v;eps_x];

V(k+1,j) = -0.5*eps_tilde.'*Pi_tilde*eps_tilde;
end
end

%% Surface plot 
trim = 20;
cut_off_range = trim:model.nt-trim;

Free_energy_curve = figure;
hold on
V_surf = surf(X_varied,model.T(cut_off_range),V(cut_off_range,:));
V_surf.EdgeAlpha = 0.1;
p_DEM = plot3(x_DEM(2,cut_off_range),model.T(cut_off_range),...
    V_real(cut_off_range),'LineWidth',2);
legend([V_surf,p_DEM],{'Free Energy Curves','DEM State Estimate'}...
    ,'interpreter','latex','Location','NorthWest')
view(45,50)
caxis(1.0e+08*[-2 -0.0001])
ax = gca;
ax.FontSize = 15;


zlabel('VFE','Interpreter','latex')
xlabel('Varied $\dot \phi$[rad]','interpreter','latex')
ylabel('Time[s]','interpreter','latex')

%% Save the figure
saveas(Free_energy_curve,'Figures/Free_energy_curve.eps','epsc2')
saveas(Free_energy_curve,'Figures/Free_energy_curve.jpg','jpg')
saveas(Free_energy_curve,'Figures/Free_energy_curve.fig','fig')