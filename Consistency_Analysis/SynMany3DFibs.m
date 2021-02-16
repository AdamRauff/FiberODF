% Adam Rauff
% MRL
% 4/20/2019

% This function generates a syhtetic image with fibers, return the image
% and its exact ODF, and then approximates the fiber ODF from the image.

function [ODF, GFA] = SynMany3DFibs(D, theta, phi, stdev, ns, lmin, lmax, TR, Nsph, Ncart, sig)

%allocate memory for image
I = zeros(D,D,D,'logical');

SampODF = zeros(size(Nsph,1),1);

% create synthetic image
[I, ~] = SeedFibers(I, ns, theta, phi, stdev, lmax, lmin, D, Ncart, SampODF, sig);

% set flags for Fiber3D function
params.saveSpec = false;
params.VisRadSpec = false;
params.VisSphPS = false;
params.VisIfft = false;

params.Thigh = 1.0;
params.Tlow = 20;

params.Ncart = Ncart;
params.Nsph = Nsph;
params.TR = TR;

% run qball to evaluate ODF from I
[ODF, GFA] = Fiber3D(I, 60, params);

% perform min-max normalization scaled by GFA



% compute generalized fractional anisotropy
% exGFA = compGFA(exODF);

%save([savedir, filesep, filename],'I', 'dist','ODF','GFA','exODF','exGFA','Nsph','Ncart');
end

