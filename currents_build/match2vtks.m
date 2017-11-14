%Match SRC shape to TGT shape. Shapes are vtk file. Saves the final result
%as well as registration evolution. OUTNAME is the output filename without
%extension. SIGMAW and SIGMAV are the kernel width of currents
%representation and diffeomorphisms (transform rigidity) respectively.
%
%   MATCH2VTKS(SRC, TGT, OUTNAME, SIGMAW, SIGMAV)

function s = match2vtks(srcsimple, tgtfilename, outname,nmins,ratio)

tic
% Opening files
[src.x src.vx] = read_vtk_shape3D(srcsimple);
[tgt.x tgt.vx] = read_vtk_shape3D(tgtfilename);


% Registration parameters
%original
% transfo.x                 = src.x';
% transfo.vx                = src.vx';
% transfo.usegrid           = 1;
% transfo.ratio             = 0.3;
% transfo.sigmaV            = sigmaV;
% transfo.stdV              = 1;
% transfo.gammaR            = 0.001;
% transfo.rigidmatching     = 0;
% transfo.numbminims        = 1;
% transfo.optim_verbosemode = 1;
% transfo.optim_maxiter     = 50;

% transfo.numbminims        = 3;
% transfo.optim_maxiter     = 1e-4;

% Registration parameters
% nico's suggestions
transfo.x                 = src.x';
transfo.vx                = src.vx';
transfo.usegrid           = 1;
transfo.ratio             = 0.3;
transfo.sigmaV            = 'auto';
transfo.stdV              = 1;
transfo.gammaR            = 0;
transfo.rigidmatching     = 0;
transfo.numbminims        = str2num(nmins);
transfo.optim_verbosemode = 1;
transfo.optim_maxiter     = 500;
transfo.optim_breakratio  = str2num(ratio);


t = cell(1,1);
t{1}.method  = 'surfcurr';
t{1}.y       = tgt.x';
t{1}.vy      = tgt.vx';
t{1}.vx      = src.vx';
t{1}.sigmaW  = 5;
t{1}.usegrid = 1;
t{1}.ratio   = 0.3;



% Do the registration
transfo = match(transfo, t);
fprintf('dist(id, phi)=%f    nNiter=%d\n', transfo.distIdPhi, length(transfo.J));

res.x  = transfo.X(:, :, transfo.T);
res.vx = src.vx';

% Set output filename prefix
f = [outname];

% Save the transformation
disp(['Saving transformation: ' f '.mat']);
save([f '.mat'], 'transfo');

% Save registered mesh
disp(['Saving registered mesh: ' f '.vtk']);
save_vtk_shape3D([f '.vtk'], res.x', res.vx', []);

% Save registration evolution
%disp('Saving registration evolution');
%for I=1:size(transfo.X,3)
%   save_vtk_shape3D([f '_' num2str(I) '.vtk'], transfo.X(:,:,I)', res.vx', []);
%end
toc
s=1;
end