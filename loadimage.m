function images=loadimage(path,filelist)
%load color image and change to cb channel, save in cell

N=length(filelist);
images=cell(N,1);
for i=1:N
    rgbimg=imread([path filesep filelist{i}]);
    rgbimg=im2uint8(rgbimg);
    ycbcr=rgb2ycbcr(rgbimg);
    cb=ycbcr(:,:,2);
    images{i}=cb;
end
