% Adam Rauff
% 4/3/2019
% MRL - angiogenesis project
% This function implements the Q-Ball algorithm developed in the 
% diffusion MRI literature.
% for relevant papers, see:
% 1) Tuch, Magnetic Resonance in Medicine, 2004
% 2) Hess et al, Magnetic Resonance in Medicine, 2006
% 3) Descoteaux, Magnetic Resonance in Medicine, 2007

function [ODF] = qball(sph_Amp, Nsph, basisOrder)

% FRT diagonal matrix with 0-order legendre polymoials evaluated at 0.
C = compLapBel_Coef(basisOrder);
 
% construct spherical harmonics matrix (see descoteaux)
T = compSH(basisOrder, Nsph(:,2), Nsph(:,1)); 
Tcond = cond(T);
disp(['Condition number of SH Matrix: ',num2str(Tcond)]); 
%[C, T] = CheckSHMat(Tcond, Nsph, C, T);
 
A = T*diag(C)*((T'*T)\T');

% load precomputed FRT transform matrix computed with basisOrder = 80
% addpath([pwd, filesep, 'Pre_Computed_SH']);
% load('SH_L80_SphPts10242.mat');

ODF = A*sph_Amp;
end

