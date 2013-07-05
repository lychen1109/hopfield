function [cbspmod]=histhacksp_run(cbtarget,cbsp,K,T,shift)
%run histhacksp to get result

N=length(cbsp);
cbspmod=cell(N,1);

parfor i=1:N
   imtarget=double(cbtarget{i});
   img=double(cbsp{i});
   selection=nodeselectionsp(img,K,shift);
   tic;
   [img]=histhacksp(img,imtarget,K,T,selection);
   tt=toc;
   fprintf('No. %d image processed in %g sec\n',i,tt);
   cbspmod{i}=uint8(img);
end


