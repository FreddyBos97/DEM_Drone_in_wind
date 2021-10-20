%% Determine the noise std of the process and measurement noise and 
% construct the histograms for the Laplace assumption 
% Figures 3 a b and c
clear all
close all
clc

%% Load the data and obtain the noise vectors and characteristics for the 
%% expeiriments with No Wind (NW) and With Wind (WW)

% select files and start and end time
start_time  = 400;
end_time    = 800;

Data_files_no_wind = [1,3,5,7];
Data_files_wind = [2,4,6,8];

for i = 1:4
data_file_no_wind{i}   = Data_files_no_wind(i);
data_file_with_wind{i} = Data_files_wind(i);

Data_no_wind{i} = load_data(data_file_no_wind{i},start_time,end_time);
Data_with_wind{i} = load_data(data_file_with_wind{i},start_time,end_time); 

Data_no_wind{i} = get_model_white_box(Data_no_wind{i},0);
Data_with_wind{i} = get_model_white_box(Data_with_wind{i},0);

Data_no_wind{i} = get_noise_charact(Data_no_wind{i},0,0);
Data_with_wind{i} = get_noise_charact(Data_with_wind{i},0,0);

% Check the standard deviations of the nosie
std_NM_w_phi(i)    = std(Data_no_wind{i}.w(1,:));
std_NW_w_phid(i)   = std(Data_no_wind{i}.w(2,:));
std_WW_w_phi(i)    = std(Data_with_wind{i}.w(1,:));
std_WW_w_phid(i)   = std(Data_with_wind{i}.w(2,:));

std_NM_phi(i)    = std(Data_no_wind{i}.x_meas(1,:));
std_NW_phid(i)   = std(Data_no_wind{i}.x_meas(2,:));
std_WW_phi(i)    = std(Data_with_wind{i}.x_meas(1,:));
std_WW_phid(i)   = std(Data_with_wind{i}.x_meas(2,:));

end

%% STD Table
[mean(std_NM_w_phi),mean(std_NW_w_phid),mean(std_WW_w_phi),mean(std_WW_w_phid)]
[mean(std_NM_phi),mean(std_NW_phid),mean(std_WW_phi),mean(std_WW_phid)]

%% Plot the state and noise vectors
Plot_colors = [0, 0.4470, 0.7410
    0.8500, 0.3250, 0.0980
    0.9290, 0.6940, 0.1250
    0.4940, 0.1840, 0.5560];

%% Autocorrelation plot for With Wind
auto_lags = 50;
font_size = 15;

auto_corr_WW = figure;

for l = 1:4    
subplot(2,1,1)
hold on

a_c_plot_nw(l,:) = autocorr(Data_with_wind{l}.w(1,:),auto_lags);
plot(linspace(0,auto_lags,auto_lags+1),a_c_plot_nw(l,:),'.-','MarkerSize',15)
plot([linspace(0,auto_lags,auto_lags+1);linspace(0,auto_lags,auto_lags+1)],...
    [zeros(1,auto_lags+1);a_c_plot_nw(l,:)],'Color',Plot_colors(l,:))


ylabel('$w_{\phi}$', 'Interpreter','latex')
title(' ')
ax = gca; 
ax.FontSize = font_size;

subplot(2,1,2)
hold on
a_c_plot_ww(l,:) = autocorr(Data_with_wind{l}.w(2,:),auto_lags);
p(l) = plot(linspace(0,auto_lags,auto_lags+1),a_c_plot_ww(l,:),'.-','MarkerSize',15);
plot([linspace(0,auto_lags,auto_lags+1);linspace(0,auto_lags,auto_lags+1)],...
    [zeros(1,auto_lags+1);a_c_plot_ww(l,:)],'Color',Plot_colors(l,:))

end
legend([p(1),p(2),p(3),p(4)],{'Wind Exp 1','Wind Exp 2','Wind Exp 3','Wind Exp 4'})

main_legend = legend([p(1),p(2),p(3),p(4)],{'Wind Exp 1','Wind Exp 2','Wind Exp 3','Wind Exp 4'});
newPosition = [0.65 0.25 0.25 0.25];
newUnits = 'normalized';
set(main_legend,'Position', newPosition,'Units', newUnits);

ylabel('$w_{\dot{\phi}}$', 'Interpreter','latex')
title(' ')
ax = gca; 
ax.FontSize = font_size;

%% Histograms for with and without wind
hist_bins = 15;
font_size = 15;
hist_NW = figure;
hold on

for j = 2
    % Remove outliers
    idx_out = abs(Data_no_wind{j}.w(1,:))>=2e-3;
    Data_no_wind{j}.w(1,idx_out) = 0;

    subplot(1,2,1)
    hold on
    p1(:,j) = histfit(Data_no_wind{j}.w(1,:),hist_bins);
    xlabel('$w_{\phi}$', 'Interpreter','latex')
    ax = gca; 
    ax.FontSize = font_size;

    subplot(1,2,2)
    hold on
    p2(:,j) = histfit(Data_no_wind{j}.w(2,:),hist_bins);
    xlabel('$w_{\dot{\phi}}$', 'Interpreter','latex')
    ax = gca; % current axes
    ax.FontSize = font_size;

end
%% 

hist_WW = figure;
hold on

for j = 2   
    subplot(1,2,1)
    hold on
    p3(:,j) = histfit(Data_with_wind{j}.w(1,:),hist_bins);
    p3(1,j).DisplayName = ['Wind Exp ',num2str(j)];
    p3(1,j).FaceColor = Plot_colors(1,:);
    p3(1,j).FaceAlpha = 1;
    p3(1,j).EdgeColor = 'k';%Plot_colors(1,:);
    p3(2,j).Color = 'r'%Plot_colors(1,:);
    xlabel('$w_{\phi}$', 'Interpreter','latex')
    ax1 = gca; 
    ax1.FontSize = font_size;
    xlim([-0.0021,0.0021])
    
    subplot(1,2,2)
    hold on
    p4(:,j) = histfit(Data_with_wind{j}.w(2,:),hist_bins);
    p4(1,j).DisplayName = ['Wind Exp ',num2str(j)];
    p4(1,j).FaceColor = Plot_colors(1,:);
    p4(1,j).FaceAlpha = 1;
    p4(1,j).EdgeColor = 'k';%Plot_colors(1,:);
    p4(2,j).Color = 'r';%Plot_colors(1,:);

    xlabel('$w_{\dot{\phi}}$', 'Interpreter','latex')
    xlim([-0.25,0.25])
    ax2 = gca; 
    ax2.FontSize = font_size;

end

%% Save figures
saveas(auto_corr_WW,'Figures/auto_corr_WW.eps','epsc')
saveas(auto_corr_WW,'Figures/auto_corr_WW.jpg','jpg')
saveas(auto_corr_WW,'Figures/auto_corr_WW.fig','fig')

% saveas(auto_corr_NW,'Figures/auto_corr_NW.jpg','jpg')
% saveas(auto_corr_WW,'Figures/auto_corr_WW.jpg','jpg')

saveas(hist_NW,'Figures/hist_WM0.eps','epsc')
saveas(hist_NW,'Figures/hist_WM0.jpg','jpg')
saveas(hist_NW,'Figures/hist_WM0.fig','fig')

saveas(hist_WW,'Figures/hist_WM2.eps','epsc')
saveas(hist_WW,'Figures/hist_WM2.jpg','jpg')
saveas(hist_WW,'Figures/hist_WM2.fig','fig')
