function selection=nodeselectionsp(spimg,K,shift)
%Select nodes which can be processed together on spatial domain
%The proper range for Cb and Cr value is [16 240]

M=size(spimg,2);
if mod(M,5)>=shift
    cols=shift:5:floor(M/5)*5+shift;
else
    cols=shift:5:(floor(M/5)-1)*5+shift;
end

selection=zeros(size(spimg));
selection(:,cols)=one;
candidate=spimg(selection==1);
inrange=ones(size(candidate));
inrange(candidate+K>240)=0;
inrange(candidate-K<16)=0;
selection(:,cols)=inrange;

