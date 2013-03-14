function selection=nodeselection(spimg)
%Select nodes which can be processed together

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
for i=1:M
    j=1;
    while j<=M
        if bdctimg(i,j)>1
            selection(i,j)=1;
            j=j+5;
        else
            j=j+1;
        end
    end
end
