function [ IMNoise ] = CreateNoisyImage(IM, SNR)
% The image is 8 bits, so the 0 < NoiseVariance < 255

%% Read in image
% number of z slices was retrieved experimentally by using imread('image
% address', indx). indx =136 works. 137 errors
% IM = zeros(168,166,136);
% 
% for i = 1:136
%     IM(:,:,i) = imread('Pre-Cell1-SRRF16-t0.tif',i);
% end

% Find variance of the noise from SNR input
NoiseVariance = std(double(IM(:)))/SNR;

% Creat Noise
noise = NoiseVariance*randn(size(IM));

% add noise to image
IMNoise = double(noise + double(IM));

% 
end

