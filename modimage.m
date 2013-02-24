function bdctimg=modimage(spimg,targetimg,T)
%function that modifies image
%input params:
%M: upbound of coeff modification

%reshape images
%turn images into bdctimg
%calculate tpm
spimg=reshape(spimg,128,128);
targetimg=reshape(targetimg,128,128);
bdctimg=blkproc(spimg,[8 8],@dct2);
bdctimg=round(bdctimg);
bdctsign=sign(bdctimg);
bdctimg=abs(bdctimg);
bdcttarget=blkproc(targetimg,[8 8],@dct2);
bdcttarget=abs(round(bdcttarget));
tpmtarget=tpm1(bdcttarget,T,1);

%initialize global V
%VG=zeros(size(spimg));

%loop all possible start point
%remove V from global V
%Vout=hopfield(bdctimg,tpmtarget,startp,M);
hopfieldnet(bdctimg,tpmtarget,T);
%update global V
%add bdctsign back to bdctimg
bdctimg=bdctimg.*bdctsign;