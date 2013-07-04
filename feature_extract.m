function feat=feature_extract(root,filelist,T)
%load feature from color images

%files=dir([root filesep '*.tif']);
N=length(filelist);
feat=zeros(2*T+1,2*T+1,N);
for i=1:N
    rgbimg=imread([root filesep filelist{i}]);
    rgbimg=im2uint8(rgbimg);
    ycbcr=rgb2ycbcr(rgbimg);
    cb=double(ycbcr(:,:,2));
    %cr=double(ycbcr(:,:,3));
    %bdctimg=blkproc(image,[8 8],@dct2);
    %bdctimg=abs(round(bdctimg));
    tm=tpm1(cb,T,2);
    feat(:,:,i)=tm;
end

feat=reshape(feat,(2*T+1)^2,N)';


