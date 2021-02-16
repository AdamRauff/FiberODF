% Adam Rauff
% 4/12/2019
% MRL - angiogenesis project

% This function computes the weighted mean of n ODFs provided following the
% methodology proposed in:
% Goh, Alvina, et al. "A Nonparamteric Riemannian Framework for Processing
% High Angular Resolution Diffusion Images (HARDI)," IEEE, 2009

% ODFs is a N x M matrix where
%   N - number of discretized points on sphere
%   M - number of ODFs
function [ meanODF, VarODF, dist ] = ODFsMeanVar(ODFs, thresh)
% define gradient descent step size
eps = 0.2;
M = size(ODFs,2);
dist = zeros(M,1);
pODFs = zeros(size(ODFs));
logP = zeros(size(ODFs));
w = (1/M)*ones(M,1);
%sqaure root transform
for i = 1:M
    pODFs(:,i) = sqrt(ODFs(:,i));
end

% intialize mean ODF to one of M ODFs
meanPODF = pODFs(:,randi(M));

prevPhi = 3;
nPhi = 2;
iter = 0;
% compute mean ODF using gradient descent
while nPhi < prevPhi
    
    iter = iter + 1;
    disp(['gradient descent iteration: ',num2str(iter)]);
    % compute logmaps from current meanPODF to every other ODF
    for j = 1:M
        logP(:,j) = LOGmap(meanPODF,pODFs(:,j));
    end
    % solve for weights (w's)
    % [V, D] = eig(logP);
    % w = V(:,1);
    
    % compute tangent vector
    phi = logP*w;
    prevPhi = nPhi;
    nPhi = norm(phi);
    disp(['tangent vector norm: ', num2str(nPhi)]);
    if nPhi < thresh
        disp('norm below threshold - mean ODF converged');
        break;
    end
    meanPODF = EXPmap(eps*phi,meanPODF);
%     for j = 1:M
%         phi = phi + w(j)*logP(:,j);
%     end
end

% compute variance
for i = 1:M
    dist(i) = acos(dot(meanPODF, pODFs(:,i)));
end
VarODF = mean(dist.^2);

% undo sqaure root transform to obtain mean ODF
meanODF = meanPODF.^2;
end

function [phi] = LOGmap(p1, p2)
    % add epsilon to denominator to ensure no division by zero
    eps = 10e-12;
    
    % if the two vectors are the same, tangent is the 0 vector
    if abs(1-dot(p1,p2)) < eps
        phi = zeros(size(p1));
    else
        phi = (p2 - dot(p1,p2)*p1)./sqrt(1-dot(p1,p2)^2) *...
        acos(dot(p1,p2));
    end
end

function [psi] = EXPmap(phi, p1)
    normPhi = sqrt(dot(phi,phi));
    eps = 10e-12;
    if normPhi < eps
        psi = p1;
    else
        psi = cos(normPhi)*p1 + sin(normPhi)*phi./normPhi;
    end
end