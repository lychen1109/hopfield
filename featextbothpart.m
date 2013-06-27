function D=featextbothpart(root,imgname,T)
%load feature from bigger part of color images

rgbimg=imread([root filesep imgname]);
image=rgb2gray(rgbimg);
image=double(image);

%compare which is bigger, red or green?
edgeimg=loadedge(root,imgname);
redmask=edgeimg(:,:,1);
greenmask=edgeimg(:,:,2);

mask=zeros(size(image));
mask(redmask==200)=1;
mask(greenmask==200)=2;

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
        if mask(i,j)==mask(i,j+1) && mask(i,j+1)==mask(i,j+2) %three pixel in the same part
            num_trans=num_trans+1;
            D(diffimg(i,j),diffimg(i,j+1))=D(diffimg(i,j),diffimg(i,j+1))+1;
        end
    end
end

D=D/num_trans;









