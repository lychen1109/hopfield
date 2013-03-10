function [spimg_tests,targetimg_tests]=samplegen(spimgs,targetimgs)
%generate 32x32 samples

spimg_tests=zeros(4,32^2);
targetimg_tests=zeros(4,32^2);
for i=1:4
    spimg=spimgs(i,:);
    spimg=reshape(spimg,128,128);
    spimg_test=spimg(1:32,1:32);
    spimg_tests(i,:)=spimg_test(:)';
    
    targetimg=targetimgs(i,:);
    targetimg=reshape(targetimg,128,128);
    targetimg_test=targetimg(1:32,1:32);
    targetimg_tests(i,:)=targetimg_test(:)';
end
