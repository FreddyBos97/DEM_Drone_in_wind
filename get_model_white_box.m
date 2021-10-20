%% Function to transform the data to the model struct which can be used by 
% the observers
function model = get_model_white_box(Data,obsv)
% Get the time
model.sam_time  = Data.dt; 
model.T         = Data.time; 
model.t_end     = max(model.T);

% These are the vectors that indicate the states for the reduced model
x_vectors = [7,10];     % phi and phi dot
y_vectors = [4];        % phi

% Get reduced state space from the white box
model.A = Data.sys_c.A(x_vectors,x_vectors);
model.B = Data.sys_c.B(x_vectors,:);
model.C = Data.sys_c.C(y_vectors,x_vectors);

% Get the data
model.x_meas = Data.x_measured(x_vectors,:);
model.y_meas = Data.y_measured(y_vectors,:);
model.v = Data.v;

% If the obervable option is selected the systems state are fully observed 
if obsv
    model.C = eye(size(model.A,1));
    model.y_meas = model.x_meas;
end

%% Preprocess input and B matrix
% Scale B and v, and subtract the mean of v 
min_v = min(min(model.v));
max_v = max(max(model.v));
B_fac = max_v-min_v;

model.v = model.v./B_fac;
model.v = model.v-mean(model.v.').';
model.B = model.B*B_fac;

model.sys_c = ss(model.A,model.B,model.C,[]);
model.sys_d = c2d(model.sys_c,model.sam_time);

% T starts at zero
model.T = model.T - min(model.T);

model.nt = length(model.T);
model.nx = size(model.A,1);
model.ny = size(model.C,1);
model.nv = size(model.v,1);
end