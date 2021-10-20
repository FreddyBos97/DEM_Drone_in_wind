function [model,brain] = generative_process(model,brain,if_dataset)
% Function that is used to determine the generalized system matrices that 
% are used in DEM, from A.A. Meera

Tt  = 0:model.sam_time:model.t_end-model.sam_time;

brain.At = kron(eye(brain.p+1),model.A);
brain.Bt = kron(eye(brain.p+1,brain.d+1),model.B);
brain.Ct = kron(eye(brain.p+1),model.C);

T = toeplitz(zeros(1,brain.p+1),[0 1 zeros(1,brain.p-1)]);
if brain.p==0; T=0;end
brain.Da = kron(T,eye(size(model.A,2)));
T = toeplitz(zeros(1,brain.d+1),[0 1 zeros(1,brain.d-1)]);
if brain.d==0; T=0;end
brain.Dv = kron(T,eye(brain.nv));
brain.D_A = brain.Da-brain.At;

brain.nx = size(model.A,1);
brain.ny = size(model.C,1);

% EDIT Precision of process noises, measurement noise, input noise
[brain.W0,brain.V0y,brain.V0v] = precision_matrix(brain.s,brain.Pw,...
                 brain.Pz,brain.sigma_v,brain.p,brain.d,brain.nv,brain.nx,brain.ny);
brain.V0 = blkdiag(brain.V0y,brain.V0v);
end