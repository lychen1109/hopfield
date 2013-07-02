function [img,tm]=histhacksp(img,imgtarget,K,T,selection)
%change only on spatial domain
%K: dymamic range of coefficients
%T: threshold of co-occurrence matrix
%if selection is not empty, Only change coeff in selection.

tpmopt=1;

%create target matrix
tmtarget=tpm1(imgtarget,T,tpmopt);
tm=tpm1(img,T,tpmopt);

if isempty(selection)
    selection=ones(size(img));
    selection(img+K>240)=0;
    selection(img-K<16)=0;
end

[SJ,SK]=find(selection);
pointsize=length(SJ);
sorted=randperm(pointsize);

for i=1:pointsize
    nodeidx=sorted(i);
    output=flaggen(img,tmtarget,SJ(nodeidx),SK(nodeidx),tm,T,K);
    if ~output.modified
        continue;
    end
    img(SJ(nodeidx),SK(nodeidx))=img(SJ(nodeidx),SK(nodeidx))+output.flag;
    tm=output.tm;
end

function output=flaggen(img,tmtarget,sj,sk,tm,T,K)
%calculate the best flag for current pixel
%dist_ori=norm(tm(:)-tmtarget(:));
dist_ori=sampledist(tm(:),tmtarget(:));
output.dist=dist_ori;
output.modified=false;
for i=-K:K
    if i==0
        continue;
    end
    out=tmmod2(img,tm,sj,sk,i,T);
    if ~out.changed
        continue;
    end
    tmnew=out.tm;
    %dist=norm(tmnew(:)-tmtarget(:));
    dist=sampledist(tmnew(:),tmtarget(:));
    if dist<dist_ori
        output.modified=true;
        dist_ori=dist;
        output.tm=out.tm;
        output.flag=i;
        output.dist=dist;
    end
end

function dist=sampledist(tm1,tm2)
%return normalized distance
%tm1=tm1(:)';
%tm2=tm2(:)';
%tm1=svmrescale(tm1,range);
%tm2=svmrescale(tm2,range);
%dist=sum(abs(tm1-tm2).*w(:));
dist=norm(tm1-tm2);



