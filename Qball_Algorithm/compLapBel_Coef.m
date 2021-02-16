% Adam Rauff
% 4/5/2019
% MRL - angiogenesis project

% Creates the diagonal entries of the Funk-Radon
% Transform 'C' 
 
% Code adapted from Iman Aganj.

function [C] = compLapBel_Coef(basisOrder)

% L is the place holder for the laplace-beltrami operator
% C is the C' from Descoteaux. The place holder for the coefficients of the
% ODF

C = zeros((basisOrder+1)*(basisOrder+2)/2, 1);
%L = C;
for k = 0:2:basisOrder
    for m = -k:k
        j = k*(k+1)/2 + m + 1;
        
        %L(j) = -k*(k+1);
        C(j) = ((-1)^(k/2))*prod(1:2:(k-1))/prod(2:2:k);
    end
end

C = (2*pi)*C;
end
