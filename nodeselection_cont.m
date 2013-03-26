function selection=nodeselection_cont(spimg,T,num)
%Select continous nodes with the num

%transform to bdct domain
bdctimg=blkproc(spimg,[8 8],@dct2);
bdctimg=abs(round(bdctimg));

%remove dc component
dcmark=ones(8,8);
dcmark(1,1)=0;
M=size(spimg,1);
dcmark=repmat(dcmark,[M/8 M/8]);
bdctimg=bdctimg.*dcmark;

selection=zeros(M,M);
if num==0
    selection(bdctimg>T)=1;
    return
else
    points=find((bdctimg>T)',num);
    selection(points)=1;
    selection=selection';
end

