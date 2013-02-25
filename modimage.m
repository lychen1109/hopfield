function bdctimgs=modimage(spimgs,targetimgs)
%function that modifies image

bdctimgs=zeros(size(spimgs));
T=10;
for i=1:size(spimgs,1)
    bdctimg=hopfieldnet(spimgs(i,:),targetimgs(i,:),T);
    bdctimgs(i,:)=bdctimg(:)';
end



