function Cb=hopfieldnet(spimg,targetimg)
%Hopfieldnet: calculate the best BDCT with hopfield network
%expected output bdctimg

%Initialization
A=500;
B=200;
C=500;
T=4;

%calculate tpm of spimg and targetimg
spimg=reshape(spimg,128,128);
targetimg=reshape(targetimg,128,128);
bdctimg=blkproc(spimg,[8 8],@dct2);
bdctsign=sign(bdctimg);
bdctimg=abs(round(bdctimg));
bdcttarget=blkproc(targetimg,[8 8],@dct2);
bdcttarget=abs(round(bdcttarget));
tpm=tpm1(bdctimg,T,2);
tpmtarget=tpm1(bdcttarget,T,2);

%calculate Cb
Cb=tpmtarget-tpm;
Cb=Cb(:);

