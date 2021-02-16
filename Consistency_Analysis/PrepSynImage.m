function [ Ip ] = PrepSynImage(I, NoiseVariance)

% convert to double, and modify fiber voxels to max intensity
Ip = double(I).*255;

% blurr image
Ip = imgaussfilt3(Ip, 8,'FilterSize', 3, 'Padding', 'symmetric');  

% strengthen signmal where fibers lie
% Ip(I) = Ip(I) + 75;

% add noise
noise = NoiseVariance*randn(size(Ip));
Ip = uint8(Ip + noise);
%Ip = imnoise(Ip, 'gaussian', 10, NoiseVariance);

end

