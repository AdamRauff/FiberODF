% Adam Rauff
% 1/14/2021
% MRL Utah
function [ flag ] = Write3DTiff_appended( I, filename, CompressionFlag )

% I - 3D image matrix to be saved as a multipage tiff
% filenmae - the directory location and desired tiff filename 
% CompressionFlag - compress image using LZW compression

% The filename must contain the path, filenmae, and extension. Example -
% [C://Desktop/MyImageFolder/DatasetX.tif]

% The compressionFlag determines weather the saved image will be
% compressed. LZW is a form of lossy compression, meaning some of the
% information in the original image will be lost. However, this makes
% visualization of large image datasets easier, and in some cases possible.

flag = false;

slices = size(I,3);

if CompressionFlag
    imwrite(I(:,:,1), filename,'Compression','LZW');

    for i = 2:slices
        imwrite(I(:,:,i), filename, 'WriteMode', 'append', 'Compression','LZW');
    end
else
    imwrite(I(:,:,1), filename,'Compression','none');

    for i = 2:slices
        imwrite(I(:,:,i), filename, 'WriteMode', 'append', 'Compression','none');
    end
end
flag = true;

end

