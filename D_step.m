function [DEM_t,DEM_x] = D_step(A,B,C,Y_embed,V0,W0,At,Bt,Da,Dv,...
                                nv,nx,ny,nt,p,d,t,sam_time,if_predict_y)
% Function that performs the state estimation step of DEM from A.A. Meera
                            
T = toeplitz(zeros(1,p+1),[0 1 zeros(1,p-1)]);
if p==0; T=0;end
Dy = kron(T,eye(ny));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic Expectation Maximization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1 = [kron(eye(p+1),C') zeros(nx*(p+1),nv*(d+1))];
k2 = (At-Da)';
k3 = [zeros(nv*(d+1),ny*(p+1)) -eye(nv*(d+1))];
k4 = kron(eye(p+1,d+1),B)';

D_A = Da-At;     

%TODO Change A3,A4,B3,B4
A1 = (k1*V0*blkdiag(kron(eye(p+1),-C),eye(nv*(d+1)))+ ...
      k2*W0*[Da-At -Bt]) + [Da zeros((p+1)*nx,nv*(d+1))];
%EDIT Input estimation as defined in paper Ajith 
A2 = ((k3*V0*blkdiag(kron(eye(p+1),-C),eye(nv*(d+1)))+ ...
      k4*W0*[Da-At -Bt])+ [zeros(nv*(d+1),(p+1)*nx) Dv]);
A3 = zeros(ny*(p+1));
%EDIT Input estimation in same way as output is estimated in DEM paper
A4 = Dv;

B1 = k1*V0;
B2 = k3*V0;
%EDIT y estimation
B3 = [Dy zeros(ny*(p+1),nv*(d+1))];
B4 = zeros(nv*(d+1),size(B1,2));

if if_predict_y
    state_sp = ss(blkdiag([A1;A2],A3,A4), [B1;B2;B3;B4], ...
        zeros(1,(nx+ny)*(p+1)+nv*(d+1)*2), zeros(1,ny*(p+1)+nv*(d+1)));
else
    state_sp = ss([A1;A2], [B1;B2], ...
        zeros(1,(nx)*(p+1)+nv*(d+1)), zeros(1,ny*(p+1)+nv*(d+1)));
end
state_sp = c2d(state_sp,sam_time,'zoh');

DEM_xx = zeros(nt,size(state_sp.C,2));
% Enter the input priors as first values for DEM estimation
DEM_xx(1,[nx*(p+1)+1:nx*(p+1)+nv]) = -Y_embed([ny*(p+1)+1:ny*(p+1)+nv],1).';

Vuu = -[k1 k2; k3 k4]*blkdiag(V0,W0)*[k1 k2; k3 k4]';
Vuy = [k1 k2; k3 k4]*blkdiag(V0,W0)*...
                [kron(eye(p+1),eye(ny)); zeros(nx*(p+1)+nv*(d+1),ny*(p+1))];
Vuv = [k1 k2; k3 k4]*blkdiag(V0,W0)*...
                [zeros(ny*(p+1),nv*(d+1)); -eye(nv*(d+1)); zeros(size(W0,1),nv*(d+1))];

% TODO Replace this Jacobian with Jacobian J in function generate_data
dfdu = [Vuu Vuy Vuv; ...
        zeros(size(Dy,1),size([Vuu Vuy Vuv],2));...
        zeros(size(Dv,1),size([Vuu Vuy Vuv],2))]  + blkdiag(Da,Dv,Dy,Dv);

% TODO Add an if-condition here to predict w and z too, similar to
% predicting y
for i = 2:nt
    if if_predict_y
        % EDIT Newton-Gauss method
        % EDIT: q what you want to predict, current state, Y_embed:
        % vertically concatenated with input prior
        q = [DEM_xx(i-1,1:nx*(p+1)+nv*(d+1)) Y_embed(1:ny*(p+1)+nv*(d+1),i-1)'];
        % EDIT f is the function, dfdu is derivative of f
        f = blkdiag([A1; A2],A3,A4)*q' + [B1;B2;B3;B4]*Y_embed(:,i-1);
        DEM_xx(i,:) = q + ((expm(dfdu*sam_time)-eye(size(f,1)))*pinv(dfdu)*f)';
    else
        % EDIT eq.(12) in paper
        DEM_xx(i,:)=(state_sp.A*DEM_xx(i-1,:)' + state_sp.B*Y_embed(:,i))';
    end
end
DEM_x = DEM_xx;
DEM_t = t;

end