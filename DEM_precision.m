%% Function that determines the DEM precision 
function [Pi_xx,Pi_vv,Pi_xv] = DEM_precision(brain)

Pi_xx = brain.Ct.'*brain.V0y*brain.Ct + brain.D_A.'*brain.W0*brain.D_A;
Pi_vv = brain.V0v + brain.Bt.'*brain.W0*brain.Bt;

Pi_xv = -brain.D_A.'*brain.W0*brain.Bt;
