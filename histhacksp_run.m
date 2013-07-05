function [cbspmod,tmnew]=histhacksp_run(cbtarget,cbsp,K,T,shift)
%run histhacksp to get result

N=length(cbsp);
cbspmod=cell(N,1);
tmnew=zeros(2*T+1,2*T+1,N);
parfor i=1:N
   imtarget=double(cbtarget{i});
   img=double(cbsp{i});
   selection=nodeselectionsp(img,K,shift);
   tic;
   [img,tm]=histhacksp(img,imtarget,K,T,selection);
   tt=toc;
   fprintf('No. %d image processed in %g sec\n',i,tt);
   cbspmod{i}=uint8(img);
   tmnew(:,:,i)=tm;
end

tmnew=reshape(tmnew,(2*T+1)^2,N)';
