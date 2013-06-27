function [featp]=feature_extract(root,T)
%load feature from color images

files=dir([root filesep '*.tif']);
N=length(files);
%feat=zeros(2*T+1,2*T+1,N);
featp=zeros(2*T+1,2*T+1,N);
for i=1:N
%     rgbimg=imread([root filesep files(i).name]);
%     image=rgb2gray(rgbimg);
%     image=double(image);
%     bdctimg=blkproc(image,[8 8],@dct2);
%     bdctimg=abs(round(bdctimg));
%     tm=tpm1(bdctimg,T,2);
%     feat(:,:,i)=tm;
    tm2=featextbothpart(root,files(i).name,T);
    featp(:,:,i)=tm2;
end

featp=reshape(featp,(2*T+1)^2,N)';


