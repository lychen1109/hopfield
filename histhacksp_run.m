function [cbspmod]=histhacksp_run(cbtarget,cbsp,K,T,shift)
%run histhacksp to get result

N=length(cbsp);
cbspmod=cell(N,1);

parfor i=1:N
   imtarget=double(cbtarget{i});
   img=double(cbsp{i});
   selection=nodeselectionsp(img,K,shift);
   tic;
   [img,distori,dist]=histhacksp(img,imtarget,K,T,selection);
   tt=toc;
   fprintf('No. %d image processed in %g sec, %g -> %g\n',i,tt,distori,dist);
   cbspmod{i}=uint8(img);
end


