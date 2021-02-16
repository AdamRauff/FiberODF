function [preODF, Nsph, Ncart ] = ReduceAmp(PS, N, VisSphAmp, params, SphDiscrete)
%% convert amplitude spectrum to spherical coordinates

% obtain evenly sampled spherical points. orientation points sampled from the unit 
% sphere using icosahedrons generated with Semechko functions from mathworks file exchange:
if SphDiscrete
    Ncart = params.Ncart;
    Nsph = params.Nsph;
    TR = params.TR;
else
    [TR, Nsph, Ncart] = getSpherPts(6);
end

%tic;
% find the nearest triangle vertex to each point (voxel location in fourier
% sprectrum).
nnVertID = nearestNeighbor(TR,N);
%toc

clear N TR
% construct unit sphere fourier spectrum by accumulating the results from
% nearest neighbor. i.e if two voxels were found to have the same vertex,
% that vertex is assigned the sum of their intensities.
preODF = accumarray(nnVertID, PS);

clear Input nnVertID
%% Convert back to cartesian coordinates 

if VisSphAmp
    figure; hold on;
    scatter3(Ncart(:,1),Ncart(:,2),Ncart(:,3), 100, preODF, '.'); colormap('gray'); 
    colorbar; caxis([min(preODF),max(preODF )]); axis image; view(3); grid on;
    xlabel('U ','FontSize', 24);
    ylabel('V','FontSize', 24);
    zlabel('W','FontSize', 24);
    view([60, 30]); hold off
    xlim([-1,1]); xticks([-1,0,1]);
    ylim([-1,1]); yticks([-1,0,1]);
    zlim([-1,1]); zticks([-1,0,1]);
    %title('Fourier Spectrum','FontSize',20);

    pause;
end

end

