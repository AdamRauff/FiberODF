function [ exODF ] = ODFContrib( exODF, theta, phi, sig, Ncart )

    % cartesian coordinates of direction
    x = cos(theta)*sin(phi);
    y = sin(theta)*sin(phi);
    z = cos(phi);
    
    % distance from each point on the sphere
    tempMat = ones(size(Ncart));
    tempMat(:,1) = tempMat(:,1)*x;
    tempMat(:,2) = tempMat(:,2)*y;
    tempMat(:,3) = tempMat(:,3)*z;
    Dist = acos(dot(tempMat,Ncart,2));
    DistNeg = acos(dot(-1*tempMat,Ncart,2));
    clear tempMat;

    % compute contribution to exact ODF by gaussian centered on [x,y,z]
    ODFcont = exp(-(Dist.^2)/(2*sig^2));    
    ODFcontNeg = exp(-(DistNeg.^2)/(2*sig^2));
    exODF = exODF + ODFcont + ODFcontNeg;
    clear ODFcont ODFcontNeg

end

