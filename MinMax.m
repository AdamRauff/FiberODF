% Adam Rauff
% 1/23/2021
% MRL, Utah

% perform a min-max normalization on an ODF, as described in Tuch 2004.

function [ ODF, GFA ] = MinMax( ODF )

GFA = compGFA(ODF);
% enforce the probability density function contraint (sum to 1) by crude 
% min-max normalization (As employed in Tuch  2004) 
ODF = ODF - (min(ODF) + (max(ODF)-min(ODF))*0.1*GFA);
ODF(ODF<0) = 0; % make sure function is non-negative everywhere
ODF = ODF./sum(ODF);
% AR Note - This ODF represents the probability of each area element on the
% unit sphere. This means sum(ODF) = 1. This is convenient for processing
% ODFs (i.e fisher rao distance). 

% compute anisotropy (see Tuch 2004)
GFA = compGFA(ODF);


end

