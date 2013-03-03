function [dist1,dist2]=distcal1(spimg,targetimg,bdctresult,T)
%calculate distance change of single image

tpmopt=1;
bdctsp=blkproc(spimg,[8 8],@dct2);
bdctsp=abs(round(bdctsp));
tpmsp=tpm1(bdctsp,T,tpmopt);

bdcttarget=blkproc(targetimg,[8 8],@dct2);
bdcttarget=abs(round(bdcttarget));
tpmtarget=tpm1(bdcttarget,T,tpmopt);

bdctresult=abs(bdctresult);
tpmresult=tpm1(bdctresult,T,tpmopt);

dist1=norm(tpmsp(:)-tpmtarget(:));
dist2=norm(tpmresult(:)-tpmtarget(:));
