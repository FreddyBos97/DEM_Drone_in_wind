%% Function that outpus the estimates for phi and phi dot, given the 
%% measurements, the model, the s value and the p value 
function [x_DEM,model,brain] = DEM_Estimate(model,brain)
[model,brain] = generative_process(model,brain,1);
model.Y_embed = generalized_process(model.y_meas.',model.prior_cause,...
    model.T,model.sam_time,model.ny,model.nv,model.p,model.d);

% DEM with unknown causes
[DEM_t,DEM_x] = D_step(model.A,model.B,model.C,model.Y_embed,brain.V0...
    ,brain.W0,brain.At,brain.Bt,brain.Da,brain.Dv,model.nv,model.nx,...
    model.ny,model.nt,model.p,model.d,model.T,model.sam_time,0);
x_DEM = DEM_x.';
end