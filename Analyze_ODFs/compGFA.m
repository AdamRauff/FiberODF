% Adam Rauff
% 4/6/2019
% MRL - angiogenesis project
% Computes the generalized fractional anisotropy 'GFA' from the ODF
%
% Code adapted from Iman Aganj.

function GFA = compGFA(ODF)

    GFA = std(ODF)/rms(ODF);
end


