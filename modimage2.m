function bdctimgs=modimage2(spimgs,targetimgs)
%function that modifies image, parfor version

bdctimgs=zeros(size(spimgs));
T=10;
M=sqrt(size(spimgs,2));
N=size(spimgs,1);
parfor i=1:N
    spimg=spimgs(i,:);
    spimg=reshape(spimg,M,M);
    targetimg=targetimgs(i,:);
    targetimg=reshape(targetimg,M,M);
    bdctimg=hopfieldnet(spimg,targetimg,T);
    bdctimgs(i,:)=bdctimg(:)';
end



