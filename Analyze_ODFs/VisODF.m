function [ out ] = VisODF( ODF, Ncart)

% in order to better visualize the ODF, we will display the ODF such that
% the sum of the point-wise probability values times their representative
% area adds up to 1. This is in contrast to their representation after
% min-max normalization, where the points represent the probability
% averaged over their areas. That is convenienvt for measuring FR distance.

% This way of visualization makes the resultant ODFs indpendent of the
% number of points used to sample the sphere.

% Assume the unit sphere is sampled such that every point has equal weighting
% (almost true with the dodecahedron sampling we use)
PtArea = 4*pi/size(Ncart, 1); % area on the unit sphere each point represents

VUC = sum(PtArea*ODF); % Volume under the curve of the probability function
ODF = ODF./VUC; % normalized probability function (pdf)

% check 
% if abs(1.0 - sum(ODF*PtArea)) < 0.001
%     disp('ODF point probabilities times surface area sums to 1');
% end

figure; hold on;
colormap('jet');
scatter3(Ncart(:,1), Ncart(:,2), Ncart(:,3), 110, ODF, '.');
xlabel('X','FontSize', 24);
ylabel('Y','FontSize', 24);
zlabel('Z','FontSize', 24);
% title('ODF','FontSize',24);
axis image;  colorbar; caxis([min(ODF), max(ODF)]);
view([60, 30]); hold off
xlim([-1,1]); xticks([-1,0,1]);
ylim([-1,1]); yticks([-1,0,1]);
zlim([-1,1]); zticks([-1,0,1]);
grid on; set(gca, 'FontSize',20); box on;
out = true;
end

