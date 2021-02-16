function [ D ] = VisLogFbrs(I)
D = max(size(I));
% flag = true;

% visualize synthetic image
mask = zeros(size(I), 'uint8');
mask(I~=0) = 255;

figure('Name','3D Rendering'); hold on;
fv = isosurface(mask,0.5);    % insert any value between 0 and 1 for the isovalue.

% isosurface likes gradients and evaluates difference in values.
% Since we have mask, a 3D matrix with object labled from 1 to TotLacNum,
% we want to display all objects and choose an isovalue that is smaller
% than 1, and bigger than 0

% note the color of the faces is given in RGB, can be changed per user
% preference. currently displaying some shade of green. This color
% preference can be added as a flag, or some user option defined earlier
patch(fv,'FaceColor','k', 'EdgeColor','none'); hold on;

% this is a shortcut for view(-37.5, 30) which defines the view point in
% polar coordinates from 0 view(azimuith, elevation)
view([45,30]);

% Reset axes to give a scale that allows to see the arrows
axis([0,D, 0,D, 0,D]);


% lighiting
camlight;
lighting gouraud; % interpolates lighting

% perspective
% camproj perspective; % parallel lines drawn from perspective. more closely correlates with human vision.
                     % displays parallel lines as converging

% label axis (note X, and Y are switched here!)
xlabel('X','FontSize', 24);
ylabel('Y','FontSize', 24);
zlabel('Z','FontSize', 24);
grid on;
hold on;
end

