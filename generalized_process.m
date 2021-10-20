function Y_embed = generalized_process(process_y,prior_cause,t,sam_time,ny,nu,p,d)
% Function to higher orders in the process output from A.A. Meera
nt = size(t,2);
E_y = [  1   0       0      0        0     0        0;
        [1 -1 0 0 0 0 0]/sam_time;
        [1 -2 1 0 0 0 0]/sam_time^2;
        [0 0 0 0 0 0 0];
        [0 0 0 0 0 0 0];
        [0 0 0 0 0 0 0];
        [0 0 0 0 0 0 0]]';

Y_embed = zeros(ny*(p+1)+nu*(d+1),nt); 
for i = 1:nt
    if i>p+1 && i<nt-p-1   % Use this section for low sampling time
        Y_embed(:,i) = [embed_Y(process_y',p+1,t(i),sam_time); ...
                       -embed_Y(prior_cause,d+1,t(i),sam_time)];
    else
        Y_embed(1:ny,i) = process_y(i,:)';
        Y_embed(ny*(p+1)+1:ny*(p+1)+nu,i) = -prior_cause(:,i);
    end
end
end