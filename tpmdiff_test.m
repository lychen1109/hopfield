function tpmdiff_test(img)
%test file of tpmdiff

img=reshape(img,128,128);
bdctimg=blkproc(img,[8 8],@dct2);
bdctimg=abs(round(bdctimg));
T=4;
tm=tpm1(bdctimg,T,1);
N=1000;%test time
result=zeros(N,1);

for i=1:N
    row=ceil(rand*128);
    col=ceil(rand*128);
    flag=ceil(rand*(2*T+1))-T-1;
    bdctnew=bdctimg;
    bdctnew(row,col)=bdctnew(row,col)+flag;
    tmnew=tpm1(bdctnew,T,1);
    result(i,1)=1-isequal(tmnew-tm,tpmdiff(bdctimg,row,col,flag));
end
fprintf('failed test number is %d\n',sum(result));
