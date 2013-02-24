function tpm=tpmdiff(array,row,col,flag,T)
%calculate tpm difference while modify single BDCT coeff

M=size(array,2);
tpm=zeros(2*T+1,2*T+1);
if col>2
    y1=threshold(array(row,col-2)-array(row,col-1),T)+T+1;
    y2=threshold(array(row,col-1)-array(row,col),T)+T+1;
    tpm(y1,y2)=tpm(y1,y2)-1;
    yp1=y1;
    yp2=threshold(array(row,col-1)-array(row,col)-flag,T)+T+1;
    tpm(yp1,yp2)=tpm(yp1,yp2)+1;
end
if col>1 && col<M
    y2=threshold(array(row,col-1)-array(row,col),T)+T+1;
    y3=threshold(array(row,col)-array(row,col+1),T)+T+1;
    tpm(y2,y3)=tpm(y2,y3)-1;
    yp2=threshold(array(row,col-1)-array(row,col)-flag,T)+T+1;
    yp3=threshold(array(row,col)+flag-array(row,col+1),T)+T+1;
    tpm(yp2,yp3)=tpm(yp2,yp3)+1;
end
if col<M-1
    y3=threshold(array(row,col)-array(row,col+1),T)+T+1;
    y4=threshold(array(row,col+1)-array(row,col+2),T)+T+1;
    tpm(y3,y4)=tpm(y3,y4)-1;
    yp3=threshold(array(row,col)+flag-array(row,col+1),T)+T+1;
    yp4=y4;
    tpm(yp3,yp4)=tpm(yp3,yp4)+1;
end
