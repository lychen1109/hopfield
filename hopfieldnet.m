function W=hopfieldnet(spimg,targetimg)
%Hopfieldnet: calculate the best BDCT with hopfield network
%expected output bdctimg
%todos: 1) filter out those illegal Vxi; 2) add input params which can
%shift colm

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
tpm=tpm1(bdctimg,T,1);
tpmtarget=tpm1(bdcttarget,T,1);

%calculate Cb
Cb=tpmtarget-tpm;
Cb=Cb(:);

%initialize W
W=zeros(2*T+1,25*128,(2*T+1)^2);
for row=1:128
    for col=1:25
        for f=1:2*T+1
            flag=f-T-1;
            if bdctimg(row,(col-1)*5+3)+flag<0
                continue;
            else
                D=tpmdiff(bdctimg,row,(col-1)*5+3,flag);
            end
            W(f,(row-1)*25+col,:)=reshape(D,1,1,(2*T+1)^2);
        end
    end
end


