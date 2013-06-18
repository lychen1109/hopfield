function D=featextpartial(root,imgname,T)
%load feature from bigger part of color images

rgbimg=imread([root filesep imgname]);
image=rgb2gray(rgbimg);
image=double(image);

%compare which is bigger, red or green?
edgeimg=loadedge(root,imgname);
numred=sum(sum(edgeimg(:,:,1)==200));
numgre=sum(sum(edgeimg(:,:,2)==200));
if numred>numgre
    mask=edgeimg(:,:,1);
else
    mask=edgeimg(:,:,2);
end

mask=double(mask);
mask=reshape(mask,size(image));
mask(mask~=200)=NaN;
mask(mask==200)=0;

image=image+mask;
diffimg=image(:,1:end-1)-image(:,2:end);
[N,M]=size(diffimg);
diffimg(diffimg>T)=T;
diffimg(diffimg<-T)=-T;
diffimg=diffimg+T+1;
msize=2*T+1;
D=zeros(msize,msize);
num_trans=0;

for i=1:N
    for j=1:M-1
        if ~isnan(diffimg(i,j)+diffimg(i,j+1))
            num_trans=num_trans+1;
            D(diffimg(i,j),diffimg(i,j+1))=D(diffimg(i,j),diffimg(i,j+1))+1;
        end
    end
end

D=D/num_trans;









