% Adam Rauff
% 8/8/2019
% MRL - angiogenesis

function [I, SampODF] = SeedFibers(I,nfibs, theta, phi, stdev, lmax, lmin, D, Ncart, SampODF, sig)

fibFam = 1;

% loop through fiber famiilies
for i = 1:fibFam
    %loop through and seed image
    for j = 1:nfibs

        %get angle from gauss distribution
        th = normrnd(theta, stdev);
        ph = normrnd(phi, stdev);
        
        % contribute to sampled distribution
        SampODF = ODFContrib(SampODF,th,ph,sig,Ncart);

        %randomly select length
        Linelen = lmin + rand(1)*(lmax-lmin);
        %length=lmax;

        %get random seed location
        loc = round(rand([1,3])*(0.93*D));
        x0 = loc(1);
        y0 = loc(2);
        z0 = loc(3);

        %get end point location
        xe = round(x0 + Linelen*cos(th)*sin(ph));
        ye = round(y0 + Linelen*sin(th)*sin(ph));
        ze = round(z0 + Linelen*cos(ph)); 

        %define 3D line equation
        lineq = @(t) t*[x0; y0; z0] + (1-t)*[xe; ye; ze];

        % step size
        s = 1/(Linelen+0.05*D);

        x = x0;
        y = y0;
        z = z0;
        for k = 1:ceil(1/s)

            if x>=1 && x<=D && y>=1 && y<=D && z>=1 && z<=D

                I(y, x, z) = true; 
                % row index = y-axis
                % column index = x-axis

            else
                break;
            end
            % end
            %update position
            newC = lineq(1-(s*k));
            x = round(newC(1));
            y = round(newC(2));
            z = round(newC(3));
        end 
    end
end

%normalize ODF
SampODF = SampODF/sum(SampODF);

%dilate image to thicken fiber lines
I = imdilate(I, strel('sphere',1)); % image is boolean at this point

% prep image by changing to 8-bit, blurring fibers, and adding noise.
% Use noise level (second argument) as 30 for control noise
I = PrepSynImage(I, 30);

% the variables I, and Iprep are kept separate, as it is much easier to
% visualize the boolean image I with isosurface.
end


