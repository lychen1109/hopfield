function feat=featuregen(imgcells,T)
%extract feature from image cells

N=size(imgcells,1);
opt=2;
feat=zeros(2*T+1,2*T+1,N);
for i=1:N
    img=double(imgcells{i});
    tm=tpm1(img,T,opt);
    feat(:,:,i)=tm;
end
feat=reshape(feat,(2*T+1)^2,N)';
