function [ODF, GFA ] = Fiber3D(image, SphHarm, params)


saveSpec = params.saveSpec; % save the power/amplitude spectrum of image
VisRadSpec = params.VisRadSpec; % view amount of info on each spherical shell in frequency spectrum
VisSphPS = params.VisSphPS; % view the spectrum (power/amplitude) after spherical "projection" (summation along radius)
VisIfft = params.VisIfft; % compute the inverse fft after filtering
%seePlane = false;  % view the plane perpendicular to each orientation
%seeNumVoxOnPlane = false; % view the number of voxels (on perpendicular plane) found for each direction
% histODF = params.histODF; % visualize histogram of ODF weights
% VisODFsphHeat = params.VisODFsphHeat; % view ODF in as a unit sphere with heat map
% VisODFradial = params.VisODFradial; % view ODF with radial position as probability value
% threshODF = params.threshODF; % discard ODF weights less then 5%
% frequency cutoffs - determined by crude trial
% These constants are hard coded 
Thigh = params.Thigh;
Tlow = params.Tlow;

if isfield(params, 'Ncart')
    SphDiscrete = true;
else
    SphDiscrete = false;
end

%thresh = 0.20;
% ask user for fiber thickness in pixels
% obtain max and min thickness of fibers

image = double(image);
IMsize = size(image);
%% construct matrices with image coordinates - for computational efficiency
[U, V, W] = size(image);
Rindex = floor(U/2 + 1);
Cindex = floor(V/2 + 1);
Zindex = floor(W/2 + 1);

% store x,y,z coordinates of every location in power spectrum.
% used for expediting computation by matrix operations, rather
% than for loops. x position changes with the column index. y position
% changes with the row index
[Um, Vm, Wm] = meshgrid(1:V, 1:U, 1:W); % AR Note 11/8 U and V where switched because meshgrid outputs matrix of "y rows and x columns"

% move center coordinate to center of spectrum, and store as signed 32bit
% integer matrix
Um = Um - Cindex;
Vm = Vm - Rindex;
Wm = Wm - Zindex;

% solve for radius of each voxel location
Rad = sqrt(Um.^2 + Vm.^2 + Wm.^2);
% Rad(Rindex, Cindex, Zindex) = 1; % ensure no radius value is 0. 
% Radius is used to obtain the inclination angle

%% Apply butterworth Filter 
% This filter darkens the edges of the images in a circular fashion, or spherical
% for 3D. This operation is crucial prior to FFT, as the FFT assumes a
% signal is continuous function from -inf to +inf. This is implemented on
% digital signals by replication the signal/image. This would cause a hard
% edge to form where the image domain ends and the replication begins. This
% edge is high frequency content that would corrupt the power spectrum. By
% darkening the edges, we are removing this edge, and instead the
% replication of the signal would be a smooth lon-frequency signal that
% will be easilt filtered out.

RadMax = min([U,V,W])/2;
% filter parameters
Do = RadMax-RadMax*0.2; % start "descend" 20% from edge
n = 10; % controls steepness of slope
% construct BW signal
H = 1./(1 + (Rad./Do).^(2*n)); % This equation comes from Gonzales textbook: "Digitial Image Processing" Third Edition page 273.

image = image.*H; % darken image edges
%% Apply Fast Fourier Transform to Image

disp('Applying FFT ...');
%2. Get fourier transform
Ispectrum = fftshift(fftn(image));
clear image;

CompSpectrum = 0;
if VisIfft
    CompSpectrum = Ispectrum; % save the complex spectrum to use with inverse transform later
end

disp('Obtaining Power Spectrum...');
PS = abs(Ispectrum).^2; % power spectrum

%4. Remove DC componenet (does not have a direction and does not constitute
%fibrillar structures)
if VisIfft
    DCcomp = CompSpectrum(Rindex, Cindex, Zindex);
end

PS(Rindex, Cindex, Zindex) = 0;

if saveSpec || VisIfft
    logSpec = log(1 + PS); % add one because log turns negative between [0,1]
    VisSpec = mat2gray(PS);
    ZPowerSliceLog = mat2gray(logSpec(:,:,Zindex));% XY Slice at middle
    ZPowerSlice = mat2gray(PS(:,:,Zindex));% XY Slice at middle
    YPowerSliceLog = mat2gray(squeeze(logSpec(Rindex,:,:))); % XZ Slice at middle
    YPowerSlice = mat2gray(squeeze(PS(Rindex,:,:)));
    logSpec = mat2gray(logSpec);
    pause;
end

% if this flag is not true, the below variables get flattened to hasten
% computations. This is the flag to visualize the inverse image after
% filtering, so I assume the user is in a debugging process. Beware, when this flag is true
% may slow down computations and be memory intensive.
if ~VisIfft
    N = [Um(:), Vm(:), Wm(:)];
    PS = PS(:);
    Rad = Rad(:);
    clear Um Vm Wm
end
 
%% Apply Ideal LP HP filters

disp('Filtering Spectrum ...');
% %filter
% Ispectrum = single(fft_filter3(Ispectrum, Tlow, Thigh, Rad));
[PS, CompSpectrum] = fft_Radial_filter3(PS, Tlow, Thigh, Rad, IMsize, VisRadSpec, VisIfft, CompSpectrum);
clear Rad

% Remake image from filtered FFT
if VisIfft

    PS = abs(CompSpectrum).^2; % obtain magnitude from complex values (power spectrum)
    logSpec = log(1 + PS); % add one so because log turns negative between [0,1]
    VisSpec = mat2gray(PS);
    ZPowerSliceLog = mat2gray(logSpec(:,:,Zindex));% XY Slice at middle
    ZPowerSlice = mat2gray(PS(:,:,Zindex));% XY Slice at middle
    YPowerSliceLog = mat2gray(squeeze(logSpec(Rindex,:,:))); % XZ Slice at middle
    YPowerSlice = mat2gray(squeeze(PS(Rindex,:,:)));
    logSpec = mat2gray(logSpec);
%     Ninv = [N(:,1) + Cindex, N(:,2) + Rindex, N(:,3) + Zindex];
    
    % put back DC component before taking ifft
    CompSpectrum(Rindex, Cindex, Zindex) = DCcomp;
    Fimage = mat2gray(abs(ifftn(CompSpectrum)));
    
    % AR Note
    % The resulting image will typically need to be saved and visualized in
    % a different software. We use Fluorender here in Utah. Also, the
    % resulting image will most likely contain a nasty "moire pattern"
    % looking checkerboard over it. This will need to be further teased out
    % to tell what is left in the image. In other words, what we filtered
    % out of the imgae with the filters. Since this is image analysis, and
    % not image processing, I will not go into details here.
end

%% Process Fourier Spectrum
% normalize rows so all vectors are converted to unit length

if ~VisIfft
    N = normr(N); 
else
    disp(['The next step of obtaining the condensed power spectrum can only take place without the inverse FFT flag']);
end

% reduce size of power spectrum
disp('Reducing Spectrum to Unit Sphere ...');
[sPS, Nsph, ~] = ReduceAmp(PS, N, VisSphPS, params, SphDiscrete);
clear Ispectrum;

%% Apply Qball Algorithm
% assign weights to fiber directions based on the plane perpendicular to
% the directrions - Funk Radon Transform
% we use the Qball algorithm from DW-MRI to compute the FRT

disp('Applying QBall ...');
addpath([pwd, filesep, 'Qball_Algorithm']); % files are in subdirectory

% apply algorithm
ODF = qball(sPS, Nsph, SphHarm);

% enforce the probability density function contraint (sum to 1) by crude 
% min-max normalization (As employed in Tuch  2004) 
[ODF, GFA] = MinMax(ODF);

end


