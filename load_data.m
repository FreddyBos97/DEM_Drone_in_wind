%% Function used to load the correct data file

function data_model = load_data(file_num,start_time,end_time)

files = dir(fullfile('Matlab_data_ww/','*.mat'));
file_names = cell(size(files));
for i = 1:length(files)
    file_names{i} = files(i).name;
end

data_model = load(['Matlab_data_ww/',file_names{file_num}]);
data_model = data_model.Data_model;

data_model.v            = data_model.v(:,start_time:end_time);
data_model.x_measured   = data_model.x_measured(:,start_time:end_time);
data_model.y_measured   = data_model.y_measured(:,start_time:end_time);
data_model.time         = data_model.time(:,start_time:end_time);
data_model.wind_vel     = data_model.wind_vel(:,start_time:end_time);
end