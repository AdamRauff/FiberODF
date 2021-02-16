% Adam Rauff
% 1/27/2020
% MRL, Utah

function [ distExact, ODFs, GFA ] = RunFiber3D_NPoints( IM, params, SHorder, exODF )

% pre-allocate variables
ODFs = zeros(size(params.Ncart,1),length(SHorder));
GFA = zeros(1, length(SHorder));
distExact = zeros(1, length(SHorder));

% In this script, the Qball itself is separated from the FFT part, as the
% FFT part is performed on the same image identicialy, but the Qball
% changes. So to save computational time, the FFT part is just done once.
sPS = DoFFT_n_Sphr_Proj(IM, params);

for i = 1:length(SHorder)
    % inform user what basis order is being processed (can be time
    % intensive)
    disp(' ');
    disp(['maxSH: ',num2str(SHorder(i))]);
    disp(' ');
    % run qball algorithm
      
    % apply algorithm
    ODFs(:,i) = qball(sPS, params.Nsph, SHorder(i));
    
    % min-max normalization
    [ODFs(:,i), GFA(i)] = MinMax(ODFs(:,i));
    
    distExact(i) = computeFisherRao(ODFs(:,i), exODF);
    
end


end

