function [ IM ] = PreAllocateImage( IMinfo , NumFiles )

if NumFiles ~= 0
    NumFiles = length(IMinfo);
end

% determine image dimensions
Rows = IMinfo(1).Height;
Cols = IMinfo(1).Width;
IMclass = IMinfo(1).BitDepth;

% determine image type
if IMclass==8
    IMclass = 'uint8';
elseif IMclass==16
    IMclass = 'uint16';
else
    disp('Could not Identify Image class');
    error('Bit depth must be 8 or 16');
end

% allocate space for image
IM = zeros(Rows, Cols, NumFiles, IMclass); 
    
end

