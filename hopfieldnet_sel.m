function [bdctimg]=hopfieldnet_sel(spimg,targetimg,T,selection,opt)
%Hopfieldnet: calculate the best BDCT with hopfield network
%intput images are square and multiple of 8
%In this version, modify all the components in selection at once

%reshape images
imgwidth=sqrt(length(spimg(:)));
spimg=reshape(spimg,imgwidth,imgwidth);
targetimg=reshape(targetimg,imgwidth,imgwidth);

%Initialization
A=opt.A;
B=opt.B;
C=opt.C;
u0=opt.u0;
M=1;%maximum modification of coeff
flags=-M:-1;
flags=cat(2,flags,1:M);

%calculate bdct and tpm of spimg and targetimg
bdctimg=blkproc(spimg,[8 8],@dct2);
bdctimg=round(bdctimg);
bdctsign=sign(bdctimg);
bdctsign(bdctsign==0)=1;
bdctimg=abs(bdctimg);

bdcttarget=blkproc(targetimg,[8 8],@dct2);
bdcttarget=abs(round(bdcttarget));
tpmtarget=tpm1(bdcttarget,T,1);
numselection=sum(sum(selection));
[compx,compy]=find(selection);  
%calculate Cb
tpm=tpm1(bdctimg,T,1);
Cb=tpmtarget-tpm;
Cb=Cb(:);
fprintf('current Cb %g\n',norm(Cb));

%initialize W
W=zeros(2*M+1,numselection,(2*T+1)^2);
for i=1:numselection
    for f=flags
        D=tpmdiff(bdctimg,compx(i),compy(i),f,T);
        W(f+M+1,i,:)=reshape(D,1,1,(2*T+1)^2);
    end
end

W=reshape(W,(2*M+1)*numselection,(2*T+1)^2);

%initalization of Tmat
% network=ones(2*M+1,numselection);
network_size=(2*M+1)*numselection;
% [nodex,nodey]=find(network);
%Tmat=zeros(network_size,network_size);

T3=W*W';
T3=sparse(T3);
T11=[0 1 1;1 0 1;1 1 0];
diagcell=cell(1);
diagcell{1}=T11;
diagcell=repmat(diagcell,numselection,1);
T1=blkdiag(diagcell{:});
T1=sparse(T1);
Tmat=-B*ones(network_size,network_size);
Tmat=-A*T1+Tmat-C*T3;

%calculate I
N=numselection;
I=N*B+C*W*Cb;

%Calculate u00, U, V, and E
u00=u0*artanh(2/(2*M+1)-1);
U=ones(network_size,1)*u00+(rand(network_size,1)*0.2-0.1)*u0;
%         U=ones(2*M+1,bwidth*L)*(-1)*u0;
%         U(M+1,:)=u0;
%         U=U(vmask);
V=nodeg(U,u0);
E=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
fprintf('iter:0  E0=%g\n',E);
%[fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode,bwidth);
%fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n',fall,f1,f2,f3);

%sequential update of nodes
iter=0;
while 1
    iter=iter+1;
    idx=randperm(network_size);
    for i=1:network_size
        m=idx(i);
        U(m)=Tmat(m,:)*V(:)+I(m);
        V(m)=nodeg(U(m),u0);
    end
    Enew=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
    fprintf('iter:%d  Enew=%g\n',iter,Enew);
    %[fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode,bwidth);
    %fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n',fall,f1,f2,f3);
    %stack=stafun(abs((Enew-E)/E),stack);
    %fprintf('length of stack %d\n',length(stack));
    %fprintf('mean of stack %g\n\n',mean(stack));
    if Enew>=E
        break;
    else
        E=Enew;
    end
end

%validation of V
V=reshape(V,2*M+1,numselection);
Vcheck=sort(V,1,'descend');
fprintf('mean of largest element of each column is %g\n',mean(Vcheck(1,:)));
fprintf('mean of second largest element of each column is %g\n',mean(Vcheck(2,:)));
V(V>0.5)=1;
V(V<=0.5)=0;

if ~isequal(sum(V),ones(1,numselection))
    fprintf('Not permutation matrix\n');
end

Sout=repmat((-M:M)',1,numselection);
Vout=sum(V.*Sout);
Vglobal=zeros(size(selection));
Vglobal(selection==1)=Vout;
newtpm=tpm1(bdctimg+Vglobal,T,1);
dist=norm(newtpm(:)-tpmtarget(:));
if dist<norm(Cb)
    fprintf('current distance %g\n',dist);
else
    fprintf('new distance is worse than default\n');
end
bdctimg=(bdctimg+Vglobal).*bdctsign;

%objective function
% function [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,network_size,nodex,nodey)
% f1=0;
% for m=1:network_size
%     for n=1:network_size
%         if nodey(m)==nodey(n) && nodex(m)~=nodex(n)
%             f1=f1+V(m)*V(n);
%         end
%     end
% end
% f1=f1*A/2;
% 
% f2=B/2*(sum(V(:))-N)^2;
% 
% f3=repmat(V,[1 size(W,2)]).*W;
% f3=sum(f3);
% f3=sum((Cb-f3').^2);
% f3=C/2*f3;
% fall=f1+f2+f3;

%stack function: add newvalue into stack
% function stack=stafun(newvalue,stack)
% l=length(stack);
% if l<3
%     stack(l+1)=newvalue;
% else
%     stack(1:2)=stack(2:3);
%     stack(3)=newvalue;
% end






    




    













