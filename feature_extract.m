function feat=feature_extract(root,T)
%load feature from color images

%root='C:\data\uncompre\4cam_auth';
files=dir([root filesep '*.tif']);
N=length(files);
feat=zeros(N,(2*T+1)^2);
for i=1:N
    rgbimg=imread([root filesep files(i).name]);
    image=rgb2gray(rgbimg);
    image=double(image);
    bdctimg=blkproc(image,[8 8],@dct2);
    bdctimg=abs(round(bdctimg));
    tm=tpm1(bdctimg,T,2);
    feat(i,:)=tm(:)';
end


