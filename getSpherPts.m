function [TR, N, Ncart] = getSpherPts(subDiv)
% orientation points were sampled for the unit sphere using icosahedrons
% generated with Semechko functions from mathworks file exchange:
% https://www.mathworks.com/matlabcentral/fileexchange/
% 37004-suite-of-functions-to-perform-uniform-sampling-of-a-sphere
%addpath([pwd,filesep,'S2 Sampling Toolbox']);

% 1) a mesh was created with IcosahedronMesh function
TR = IcosahedronMesh;

% 2) Mesh was subdivided 6 times to create a fine mesh 
TR = SubdivideSphericalMesh(TR,subDiv);

Ncart = TR.Points;
% N = Pts((Pts(:,3)>=0),:);
[az, el, ~] = cart2sph(Ncart(:,1), Ncart(:,2), Ncart(:,3));

N = [az, el]; clear az el;

% 4) radians converted to degrees
% N(:,1) = rad2deg(N(:,1));
% N(:,2) = rad2deg(N(:,2));

N(:,2) = pi/2 - N(:,2); % convert elevation to inclination from +z axis
%AmpSph = N;
% 5) moved range of azimuth to [0,360]


end

