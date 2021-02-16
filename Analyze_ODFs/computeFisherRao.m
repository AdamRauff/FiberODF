% Adam Rauff
% MRL - angiogenesis
% 4/11/2019
% This function computes the distance between 2 Orientation Distribution
% Functions (ODFs) using the Fisher-Rao metric.

% Goh, Alvina, et al. "A Nonparametric Riemannian Framework for Processing
% High Angular Resolution Diffusion Images (HARDI)" IEEE, 2009

% this function return the distance in radiancs between 2 ODFs
function [ distDeg ] = computeFisherRao(ODF1, ODF2)

% The formulation for this assumes that functions are continuous unit
% vector that sum to 1. This the leads the square root transform to yield
% unit vectors over the unit Hilbert sphere. Thus we need to multiply our
% ODFs by a constant so that the point-wise values represent the functional
% values over the area of the sphere they represent
% ODF1 = ODF1*PtArea;
% ODF2 = ODF2*PtArea;
% 
% % Crude min-max normalization, As employed in Tuch et al 2004
% ODF1 = ODF1-min(ODF1);
% ODF1 = ODF1/sum(ODF1);
% 
% ODF2 = ODF2-min(ODF2);
% ODF2 = ODF2/sum(ODF2);

% sqaure root transform
p1 = sqrt(ODF1);
p2 = sqrt(ODF2);


% compute tanget vectors
% distance
distRad = acos( dot(p1,p2) );
distDeg = rad2deg(distRad);
end

