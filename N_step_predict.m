function x_predict = N_step_predict(model,v_est,N_step,T_range)

NT = length(T_range);
x_meas = model.x_meas(:,T_range);

x_predict = zeros(model.nx,NT);
x_predict(:,[1:N_step-1]) = x_meas(:,1)*ones(1,N_step-1);

for j = 1:NT-N_step+1
    clear x_predict_N v_N
    x_predict_N(:,1) = x_meas(:,j);
    v_N = v_est(:,j:j+N_step-1);
    for k = 1:N_step-1
        x_predict_N(:,k+1) = model.sys_d.A*x_predict_N(:,k) +...
            model.sys_d.B*v_N(:,k);
    end
    x_predict(:,j+N_step-1) = x_predict_N(:,end);
end

end