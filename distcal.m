function [dist1,dist2]=distcal(imgoris,targetimgs,bdctresults)
%calculate result of distance change

T=10;
L=size(imgoris,1);
dist1=zeros(L,1);
dist2=zeros(L,1);
for i=1:L
    imgori=imgoris(i,:);
    imgori=reshape(imgori,128,128);
    bdctori=blkproc(imgori,[8 8],@dct2);
    bdctori=abs(round(bdctori));
    tpmori=tpm1(bdctori,T,2);
    
    targetimg=targetimgs(i,:);
    targetimg=reshape(targetimg,128,128);
    bdcttarget=blkproc(targetimg,[8 8],@dct2);
    bdcttarget=abs(round(bdcttarget));
    tpmtarget=tpm1(bdcttarget,T,2);
    
    bdctresult=bdctresults(i,:);
    bdctresult=reshape(bdctresult,128,128);
    bdctresult=abs(bdctresult);
    tpmresult=tpm1(bdctresult,T,2);
    
    dist1(i)=norm(tpmori(:)-tpmtarget(:))^2;
    dist2(i)=norm(tpmresult(:)-tpmtarget(:))^2;
end
