%this function takes in an FFT matrix, low and high period cuttoff
%values and the center coordinates of the fft matrix and returns a matrix
%with all values outside of this removed

function [PS, compSpectrum] = fft_Radial_filter3(PS, Tlow, Thigh, Rad, IMsize, VisRadSpec, VisIfft, compSpectrum)

% convert spatial period (in pixles) to frequency cuttoff 
N = min(IMsize)/2; % this is only apropriate for cubic images. for non-cubic images, this would not be a constant - filter would be ellipsoidal instead of spherical
flow = 1/Tlow;
fhigh = 1/Thigh;

F = Rad./N;
% F = N./Rad;
% find indices above and below acceptable range of frequencies
ZilchInd = (F<flow | F>fhigh);

% zero values in power spectrum
PS(ZilchInd) = 0;
%Cartcoord(ZilchInd,:) = [];

if VisIfft
    %compSpectrum = compSpectrum(:); % flatten complex power spectrum
    % zero values in power spectrum
    compSpectrum(ZilchInd) = 0;
end

% view amount of info on each spherical shell in frequency spectrum
if VisRadSpec
    % AR 4/3/2019
    % This information could be used to determine frquency cutoffs that
    % depend on the data rather then hard coded numbers.
    %Rad(ZilchInd) = [];
    maxRad = floor(fhigh * N);
    minRad = ceil(flow*N);
    RadTicks = (minRad:maxRad)';
    subs = nearestpoint(Rad(:), RadTicks);
    RadVis = accumarray(subs, PS(:));
    clear subs maxRad;
    
    figure; 
    plot(RadTicks, RadVis, 'k','LineWidth',2);
    grid on; box on;
    xlabel('Frequency [voxels]','FontSize',20);
    ylabel('Magnitude','FontSize',20);
    title('Power Spectrum Intensity over Spherical Shells','FontSize',24);
    pause;
end


end