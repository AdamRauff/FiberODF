% Adam Rauff
% 1/23/2021
% MRL

% This script test the ability of our algorithm to approximate highly
% anisotropic ODFs. That is, ODFs with high angular resolution. This is
% tested by trying to approximate the dirac detla function, which
% corresponds to parallel fibers.
% prep space
clear; clc; close all

%% add paths to needed functions
addpath([pwd, filesep, '..']);
addpath([pwd, filesep, '..', filesep, 'Qball_algorithm']);
addpath([pwd, filesep, '..', filesep, 'Analyze_ODFs']);
addpath([pwd, filesep, '..', filesep, 'S2 Sampling Toolbox']);

% exact ODF and syntehtic fibers pointing along x-axis
load('Apprx_Dir_Delt_Workspace.mat');

IM = I_8bit; clear I_8bit;

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
params.Thigh = 1.0; % ideal high pass filter cutoff (measured in pixels)
params.Tlow = 50; % ideal low pass filter cutoff (measured in pixels)
params.Ncart = Ncart_SD6;
params.Nsph = Nsph_SD6;
params.TR = TR_SD6;

maxSH_SD6 = 84;
maxSH_SD5 = 78;
maxSH_SD4 = 46;

increment = 4;

% pre-allocate variables
% spherical harmonic order L
SHorder_SD6 = 4:increment:maxSH_SD6; 
SHorder_SD5 = [4:increment:(maxSH_SD5 - 2), maxSH_SD5];
SHorder_SD4 = [4:increment:(maxSH_SD4 - 2), maxSH_SD4];

%% SD6 - 40k points
disp(' ');
disp('Beginning Analysis with 40k points');
disp(' ');
[distExact_SD6, ODFs_SD6, GFA_SD6] = RunFiber3D_NPoints(IM, params, SHorder_SD6, exODF_SD6);

%% SD5 - 10k points
% Set parameters to those of subdivision 5 - 10k points
params.Ncart = Ncart_SD5;
params.Nsph = Nsph_SD5;
params.TR = TR_SD5;

disp(' ');
disp('Beginning Analysis with 10k points');
disp(' ');
[distExact_SD5, ODFs_SD5, GFA_SD5] = RunFiber3D_NPoints(IM, params, SHorder_SD5, exODF_SD5);

%% SD4 - 2k points
% Set parameters to those of subdivision 4 - 2k points
params.Ncart = Ncart_SD4;
params.Nsph = Nsph_SD4;
params.TR = TR_SD4;
[distExact_SD4, ODFs_SD4, GFA_SD4] = RunFiber3D_NPoints(IM, params, SHorder_SD4, exODF_SD4);

%% Make figure
figure;
yyaxis left;
plot(SHorder_SD6, distExact_SD6,'k','LineWidth',2); hold on;
plot(SHorder_SD5, distExact_SD5,'k--','LineWidth',2); hold on;
plot(SHorder_SD4, distExact_SD4 ,'k-.','LineWidth',2); hold on;

yyaxis right;
plot(SHorder_SD6, GFA_SD6, 'b', 'LineWidth', 2); hold on
plot(SHorder_SD5, GFA_SD5, 'b--', 'LineWidth', 2);
plot(SHorder_SD4, GFA_SD4, 'b-.', 'LineWidth', 2);
box on; grid on; hold on;
set(gca,'FontSize',20);
