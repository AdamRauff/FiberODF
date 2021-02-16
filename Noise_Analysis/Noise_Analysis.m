% Adam Rauff
% 4/27/2019
% MRL

% this script tests for the maximum SH term used for Qball approximation of
% ODFs

% prep space
clear; clc; close all

% add paths to needed functions
addpath([pwd, filesep, '..']);
%addpath([pwd, filesep, '..', filesep, '..', filesep, 'VerifyQball', filesep, 'ScottWood']);
addpath([pwd, filesep, '..', filesep, 'Qball_algorithm']);
addpath([pwd, filesep, '..', filesep, 'Analyze_ODFs']);
addpath([pwd, filesep, '..', filesep, 'S2 Sampling Toolbox']);

%% Read in image of parallel fibers
% Synthetic images
NoiseLess = 'X_aligned_350_fibs_len_120-200_10k_pts.mat';
load(NoiseLess);

% Actin Filaments on Chondrocyte - provided by Dr. Scott Wood
SWIM_Raw = zeros(168, 166, 136);
SW_IM = zeros(136, 136, 136);

Relative_Dir = ['..', filesep,'Biomedical_Image_Data', filesep, 'SW_Actin_Chondrocyte'];
filename = 'Pre-Cell1-SRRF16-t0.tif';
for i = 1:136
    SWIM_Raw(:,:,i) = imread([Relative_Dir,filesep, filename],i);
    SW_IM(:,:,i) = SWIM_Raw(33:end, 31:end, i);
end
clear SWIM_Raw

%% Set SNR and harmonic parameters
% Signal to noise ratio
SNR = [0.2, 0.4, 0.6, 0.8, 1, 2:2:20];

n = length(SNR);
maxSNR = max(SNR);
minSNR = min(SNR);

% spherical harmonic order
SHorder = [10, 20, 30, 40, 50, 60, 70, 78]; 

%% Analyze Simulated Images

% pre-allocate matrices for memory efficiency
ODF_Sim_fibs = zeros(size(Ncart,1),1);
distNoiseLess_simFibs = zeros(length(SHorder),length(SNR));
ODF_Noiseless_simFibs = zeros(size(Ncart,1),length(SHorder));
GFA_NoiseFree = zeros(length(SHorder),1);
GFA_Noiseless_Diff = zeros(length(SHorder),length(SNR));

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
params.Thigh = 1.5; % ideal high pass filter cutoff (measured in pixels)
params.Tlow = 20; % ideal low pass filter cutoff (measured in pixels)
params.Ncart = Ncart;
params.Nsph = Nsph;
params.TR = TR;

% noiseless approximation
% FFT part of algorithm
disp('Noiseless approximations');
disp(' ');
sPS = DoFFT_n_Sphr_Proj(I_16bit, params);

% run qball to evaluate ODF from Image
for i = 1:length(SHorder)
    
    % inform user what basis order is being processed (can be time
    % intensive)
    disp(' ');
    disp(['maxSH: ',num2str(SHorder(i))]);
    disp(' ');
    % run qball algorithm
      
    % apply algorithm
     ODF_Noiseless_simFibs(:,i) = qball(sPS, params.Nsph, SHorder(i));
    
    % min-max normalization
    [ODF_Noiseless_simFibs(:,i), ~] = MinMax(ODF_Noiseless_simFibs(:,i));
end

clear ODF_temp GFA_temp

for j = 1:length(SNR)
    
    disp(' ');
        disp(['SNR: ',num2str(SNR(j))]);
    disp(' ');

    IM_Noise = CreateNoisyImage(I_16bit, SNR(j));
    sPS = DoFFT_n_Sphr_Proj(IM_Noise, params);
    
    for i = 1:length(SHorder)
        
        disp(' ');
        disp(['maxSH: ',num2str(SHorder(i))]);
        disp(' ');

         % apply algorithm
         ODF_Sim_fibs = qball(sPS, params.Nsph, SHorder(i));

        % min-max normalization
        [ODF_Sim_fibs, ~] = MinMax(ODF_Sim_fibs);
       
        % measure distance from noiseless approximation
        distNoiseLess_simFibs(i,j) = computeFisherRao(ODF_Noiseless_simFibs(:,i),ODF_Sim_fibs);
        %GFA_Noiseless_Diff(i,j) = abs(GFA_NoiseFree(i) - GFA_temp);

    end
end

clear IM_Noise
figure;

for i = 1:length(SHorder)
    plot(SNR, distNoiseLess_simFibs(i,:), 'LineWidth', 2, 'DisplayName',['L = ',num2str(SHorder(i))]); hold on
end

title('Synthetic Data','FontSize',24);
xlabel('SNR'); ylabel('Distance');
set(gca,'FontSize',18)
box on; grid on;

% Plot GFA
% figure;
% plot(SNR, GFA_Noiseless_Diff(1,:), 'k-o', 'LineWidth', 2); hold on;
% plot(SNR, GFA_Noiseless_Diff(2,:), 'r-o', 'LineWidth', 2); hold on;
% plot(SNR, GFA_Noiseless_Diff(3,:), 'b-o', 'LineWidth', 2); hold on
% plot(SNR, GFA_Noiseless_Diff(4,:), 'g-o', 'LineWidth', 2); hold on
% plot(SNR, GFA_Noiseless_Diff(5,:), 'c-o', 'LineWidth', 2); hold on
% legend('L = 10', 'L = 25', 'L = 40', 'L = 55', 'L = 70');
% xlabel('SNR'); ylabel('GFA difference ');
% set(gca,'FontSize',18)
% box on; grid on;
%% Analyze biomedical Images

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

% pre-allocate matrices for memory efficiency
ODF_SW = zeros(size(Ncart,1),1);

%GFA_SW = zeros(1,length(SNR));
distNoiseLess_SW = zeros(length(SHorder),length(SNR));
ODF_Noiseless_SW = zeros(size(Ncart,1),length(SHorder));
GFA_NoiseFree = zeros(length(SHorder),1);
GFA_Noiseless_Diff = zeros(length(SHorder),length(SNR));

disp('Analyzing Scott Wood Data'); disp('');
% noiseless approximation
% run qball to evaluate ODF from I

sPS = DoFFT_n_Sphr_Proj(SW_IM, params);

for i = 1:length(SHorder)
    
    % inform user what basis order is being processed (can be time
    % intensive)
    disp(' ');
    disp(['maxSH: ',num2str(SHorder(i))]);
    disp(' ');
    % run qball algorithm
      
    % apply algorithm
     ODF_Noiseless_SW(:,i) = qball(sPS, params.Nsph, SHorder(i));
    
    % min-max normalization
    [ODF_Noiseless_SW(:,i), ~] = MinMax(ODF_Noiseless_SW(:,i));
    
end

for j = 1:length(SNR)
    
    SW_NoiseyIM = CreateNoisyImage(SW_IM, SNR(j));
    sPS = DoFFT_n_Sphr_Proj(SW_NoiseyIM, params);
    for i = 1:length(SHorder)
        
        
        disp(' ');
        disp(['maxSH: ',num2str(SHorder(i))]);
        disp(' ');

         % apply algorithm
        ODF_SW = qball(sPS, params.Nsph, SHorder(i));

        % min-max normalization
        [ODF_SW, ~] = MinMax(ODF_SW);
       
        % measure distance from noiseless approximation
        distNoiseLess_SW(i,j) = computeFisherRao(ODF_Noiseless_SW(:,i),ODF_SW);
        
    end
end

clear SW_NoiseyIM;
figure;

for i = 1:length(SHorder)
    plot(SNR, distNoiseLess_SW(i,:), 'LineWidth', 2, 'DisplayName',['L = ',num2str(SHorder(i))]); hold on
end

title('SW Data','FontSize',24);
xlabel('SNR'); ylabel('Distance');
set(gca,'FontSize',18)
box on; grid on;

