function feat=featuregen(tmlist)
%extract feature from unnormalized transition matrix

N=size(tmlist,1);
feat=zeros(size(tmlist));
for i=1:N
    tm=tmlist(i,:);
    tm=tm/sum(tm);
    feat(i,:)=tm;
end
