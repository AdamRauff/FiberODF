function [ IM ] = ReadMultiPageTiff( pathToImage)

IMinfo = imfinfo(pathToImage);

% pre-allocate space for image
IM = PreAllocateImage(IMinfo, 0); 

% populate image with the slices that are read from file
TifLink = Tiff(pathToImage, 'r');
for i = 1:length(IMinfo)
    %disp(['Image - ',num2str(i)]);
    TifLink.setDirectory(i);
    IM(:,:,i) = TifLink.read();
end
TifLink.close();

clear I_temp;

end

