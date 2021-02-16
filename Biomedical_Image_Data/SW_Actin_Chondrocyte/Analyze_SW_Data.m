% Adam Rauff
% 1/14/2021
% MRL

% This script evaluates the effect of the number of spherical harmonic
% terms that are chosen (by the user) to represent the condensed power
% spectrum, and consequently the resulting orientation distribution (ODF).

% prep workspace
clear; clc; close all

%% add paths to needed functions
addpath([pwd, filesep, '..']);
addpath([pwd, filesep, '..', filesep, '..']);
addpath([pwd, filesep, '..', filesep, '..', filesep, 'Qball_algorithm']);
addpath([pwd, filesep, '..', filesep, '..', filesep, 'Analyze_ODFs']);
addpath([pwd, filesep, '..', filesep, '..', filesep, 'S2 Sampling Toolbox']);

%% Read in image
% number of z slices was retrieved experimentally by using imread('image
% address', indx). indx =136 works. 137 errors
IM = zeros(168,166,136);

for i = 1:136
    IM(:,:,i) = imread('Pre-Cell1-SRRF16-t0.tif',i);
end

% Choose cubical segment of image
IM = IM(33:end, 31:end, 1:136);

%% run our algorithm

% obtain sphere discretization points
[TR, Nsph, Ncart] = getSpherPts(6);


% set flags for Fiber3D function. Used for debugging. Warning, this will
% slow down the process
params.saveSpec = false;
params.VisRadSpec = false;
params.VisIfft = false;
params.VisSphPS = false;

% set parameters for Fiber3D function.
% The ideal High/low pass filters are encoded in pixels. This way the
% analyst can approximate the sturctures they would like to exlcude
% intuitively. Ex. If my fibers in image X are 4-10 pixels thick, the high
% filter cutoff should be smaller than 4 and the low filter cutoff should
% be higher than 10. It is worth noting the ringing effect can be strong
% with discrete signals, so these should only be used as ballpark
% estimates.
params.Thigh = 2.5; % ideal high pass filter cutoff (measured in pixels)
params.Tlow = 20; % ideal low pass filter cutoff (measured in pixels)
params.Ncart = Ncart;
params.Nsph = Nsph;
params.TR = TR;

maxSH = 16; % maximum number of spherical harmonic terms. 

% approximate ODF with maximum number of spherical harmonic terms
[ODF_maxSH, GFA_maxSH] = Fiber3D(IM, maxSH, params);


% pre-allocate variables
increment = 2;
SHorder = 2:increment:(maxSH - increment); % spherical harmonic order
ODFs = zeros(size(Ncart,1), length(SHorder), 'double');
GFA = zeros(1, length(SHorder));
GFA_diff = zeros(1, length(SHorder) - 1); % difference in anisotropy
distExact = zeros(1, length(SHorder)); % FR distance from "exact" ODF approximation
distPrev = zeros(1, length(SHorder) - 1);


for i = 1:length(SHorder)
    % inform user what basis order is being processed (can be time
    % intensive)
    disp(' ');
    disp(['maxSH: ',num2str(SHorder(i))]);
    disp(' ');
    % run algorithm
    
    [ODFs(:,i), GFA(i)] = Fiber3D(IM, SHorder(i), params);
    distExact(i) = computeFisherRao(ODFs(:,i), ODF_maxSH);
    
    if i > 1
        distPrev(i-1) = computeFisherRao(ODFs(:,i), ODFs(:,i-1));
        GFA_diff(i-1) = abs(GFA(i) - GFA(i-1));
    end
end

modSHorder = SHorder(2:end);
figure;
plot(SHorder,distExact, 'k', 'LineWidth', 3); hold on;
plot(modSHorder,distPrev, 'k--', 'LineWidth', 3);

yyaxis right;
plot(modSHorder, GFA_diff, 'b', 'LineWidth', 2);
legend(['Distance from L = ',num2str(maxSH)],'Distance from (L - 2)', 'Difference of Anisotropy from (L-2)');
box on; grid on; hold on;
set(gca,'FontSize',20);

