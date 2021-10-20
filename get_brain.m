%% Function to input all the system values into brain for the DEM algorithm
function brain = get_brain(model)
% Save everything into brain and model for the functions
brain.p = model.p;
brain.d = model.d;
brain.nx = model.nx;
brain.ny = model.ny;
brain.nt = model.nt;
brain.nv = model.nv;
brain.s = model.s;
brain.Pw = model.Pw;
brain.Pz = model.Pz;
brain.sigma_v = model.sigma_v;
end

