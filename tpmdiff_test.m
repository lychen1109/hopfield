function tpmdiff_test(img)
%test file of tpmdiff

img=reshape(img,128,128);
bdctimg=blkproc(img,[8 8],@dct2);
bdctimg=abs(round(bdctimg));
T=4;
tm=tpm1(bdctimg,T,1);
N=400;%test time
result=zeros(N,1);

for i=1:N
    bdctnew=bdctimg;
    tmsig=tm;
    for row=1:128
        for coli=1:25
            %col=(coli-1)*5+3;
            col=coli;
            flag=ceil(rand*(2*T+1))-T-1;
            tmsig=tmsig+tpmdiff(bdctnew,row,col,flag);
            bdctnew(row,col)=bdctnew(row,col)+flag;            
        end
    end
    tmnew=tpm1(bdctnew,T,1);
    resulti=isequal(tmnew,tmsig);
    fprintf('test result is %d\n',resulti);
    result(i,1)=1-resulti;
end
fprintf('failed test number is %d\n',sum(result));
