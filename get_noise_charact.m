%% Function to determine noise characteristics
function model = get_noise_charact(model,ms_num,optimize_on)
% Get the process noise 
w = zeros(model.nx,model.nt);
for i = 1:model.nt-1
    w(:,i+1) = model.x_meas(:,i+1) - model.sys_d.A*model.x_meas(:,i) - ...
        model.sys_d.B*model.v(:,i);
    
end
w1 = w(1,:);
w2 = w(2,:);

% Determine the noise characteristics
cov_w = cov(w1,w2);
sigma_w1 = sqrt(cov_w(1,1));
sigma_w2 = sqrt(cov_w(2,2));

%% Matlab multistart for optimisation of s, based on the SSE between the 
% autocorrelation of the noise vs the autocorrelation of a gaussian fit
if optimize_on ==1
    opts = optimoptions(@fmincon,'Algorithm','sqp');
    ms = MultiStart;

    obj_fun_w1 = @(x) fit_obj(x,model.T,w1,sigma_w1);
    problem_w1 = createOptimProblem('fmincon','objective',...
        obj_fun_w1,'x0',0.1,'lb',0,'ub',0.5,'options',opts);

    obj_fun_w2 = @(x) fit_obj(x,model.T,w2,sigma_w2);
    problem_w2 = createOptimProblem('fmincon','objective',...
        obj_fun_w2,'x0',0.1,'lb',0,'ub',0.5,'options',opts);

    [s_opt1,s_opt_sse1] = run(ms,problem_w1,ms_num);
    [s_opt2,s_opt_sse2] = run(ms,problem_w2,ms_num);

    model.s1 = s_opt1;
    model.s2 = s_opt2;
end 

%% Fill in the proper noise charactaristics 
Pw = inv(cov_w);
Pw_diag = zeros(model.nx);
for i = 1:model.nx
    Pw_diag(i,i) = 1/cov_w(i,i);
end 
model.Pw = Pw_diag;
model.w = w;

function SSE_fit = fit_obj(x,T,w,sigma_w)
    w_gauss = convolute_gaussian([x,sigma_w],T);
    auto_corr = autocorr(w_gauss,100);
    auto_w = autocorr(w,100);
    SSE_fit = determine_sse(auto_corr,auto_w,0);
end
end 