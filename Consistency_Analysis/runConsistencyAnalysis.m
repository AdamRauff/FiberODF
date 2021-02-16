% Adam Rauff 
% 4/6/2019
% MRL - angiogenesis project

% prep space
clear; clc; close all

% add paths to needed functions
addpath([pwd, filesep, '..']);
addpath([pwd, filesep, '..', filesep, 'Qball_algorithm']);
addpath([pwd, filesep, '..', filesep, 'Analyze_ODFs']);
addpath([pwd, filesep, '..', filesep, 'S2 Sampling Toolbox']);

%specify image dimensions and name
D = 400;

% specify number of monte-Carlo runs
n = 20;

% obtain sphere discretization points to construct exact ODF
[TR, Nsph, Ncart] = getSpherPts(6);

NumSphPts = size(Ncart,1);
%ODFs = zeros(size(Ncart,1),n);
%exODFs = zeros(size(Ncart,1),n);

%specify center angle and st. deviation (in radians)
theta = 0;
phi = pi/2; 
stdEv = [0, 1, 2, 3, 5, 8];
stdrad = deg2rad(stdEv);

% ODF blur
sig = deg2rad(2); 

%specify number of seeds and max length
ns = 300;
lmax = 0.3*D;
lmin = 0.1*D;

% exact ODF
% exODF = zeros(size(Nsph,1),1);
% exODF = ODFContrib(exODF, theta, phi, stdrad, Ncart);
% exODF = exODF/sum(exODF);

%dist = zeros(n,1);
ODFs = zeros(NumSphPts,n);
%SeampODFs = zeros(NumSphPts,n);
GFA = zeros(n,length(stdrad));

thresh = 10e-8; % threshold for convergence
Distances_Degrees = zeros(n,length(stdrad));
for i = 1:length(stdrad)
    Fiber_Deviation = stdrad(i);
    for k = 1:n

        disp(['Iteration: ', num2str(k)]);

        [ODFs(:,k), GFA(k,i)] = SynMany3DFibs(D, theta, phi, Fiber_Deviation, ns, lmin, lmax, TR, Nsph, Ncart, sig); 
    end


    % compute mean ODF
    
    [mODF, VarODF, Distances] = ODFsMeanVar(ODFs, thresh);

    Distances_Degrees(:,i) = rad2deg(Distances);
    stdODF = rad2deg(sqrt(VarODF));

end