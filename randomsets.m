function featnew=randomsets(featau,featcasia)
%get random samples from combined sets

N=921;%random select 921 authentic images
idx=randperm(921+800);
featnew=[featau;featcasia];
featnew=featnew(idx(1:921),:);
selection=idx(1:921);
numau=sum(selection<=921);
numcasia=sum(selection>921);
fprintf('number of original sets %d\n',numau);
fprintf('number of new sets %d\n',numcasia);
