% Adam Rauff
% 1/26/2020

% This script is written to analyze data sets of interest are from the Focused Ion Beam
% Scanning electron microscopy (FIB-SEM) from the study of Sean Reese. The code begins by
% resampling images with anisotropic voxel aspect to isotropic voxels of
% 6 x 6 x 6 [um]. Mostly that aspect ratio. Sometimes something else

% prepare workspace
clc; clear; close all;

%% add paths to needed functions
addpath([pwd, filesep, '..']);
addpath([pwd, filesep, '..', filesep, '..']);
addpath([pwd, filesep, '..', filesep, '..', filesep, 'Qball_algorithm']);
addpath([pwd, filesep, '..', filesep, '..', filesep, 'Analyze_ODFs']);
addpath([pwd, filesep, '..', filesep, '..', filesep, 'S2 Sampling Toolbox']);


%% compose image
Block_1_FileName = 'IM_Block_Analyzed_1.tif';
Block_2_FileName = 'IM_Block_Analyzed_2.tif';

% pre-allocate space for image
IMBlock1 = ReadMultiPageTiff([pwd, filesep, Block_1_FileName]);  
IMBlock2 = ReadMultiPageTiff([pwd, filesep, Block_2_FileName]);

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
% estimates. cahrefull analysis should checkout the IFFT flag and inspect
% the power spectrum visually
params.Thigh = 2.0; % ideal high pass filter cutoff (measured in pixels)
params.Tlow = 30; % ideal low pass filter cutoff (measured in pixels)
params.Ncart = Ncart;
params.Nsph = Nsph;
params.TR = TR;

maxSH = 16; % maximum number of spherical harmonic terms. 

%% Image Block 1
% run our algorithm
[ODF, GFA] = Fiber3D(IMBlock1, maxSH, params);

%% Image Block 2

% run our algorithm
[ODF2, GFA2] = Fiber3D(IMBlock2, maxSH, params);

dist = computeFisherRao(ODF, ODF2);
