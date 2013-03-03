function bdctimg=hopfieldnet(spimg,targetimg,T)
%Hopfieldnet: calculate the best BDCT with hopfield network
%intput images are square and multiple of 8
%todo: 1. change program so that it can process image wth different size
%      2. record feature distance of every iteration

%Initialization
A=1500;
B=400;
C=20;
u0=300;
M=3;%maximum modification of coeff
L=4;%Process L lines a time

stack=[];%used in calculation of average delta Energy function
tol=5e-4;%tolerance of minimum delta Energy function

%calculate bdct and tpm of spimg and targetimg
bdctimg=blkproc(spimg,[8 8],@dct2);
bdctimg=round(bdctimg);
bdctsign=sign(bdctimg);
bdctsign(bdctsign==0)=1;
bdctimg=abs(bdctimg);
Vglobal=zeros(size(bdctimg));

bdcttarget=blkproc(targetimg,[8 8],@dct2);
bdcttarget=abs(round(bdcttarget));
tpmtarget=tpm1(bdcttarget,T,1);

[MI,NI]=size(bdctimg);
RV=MI/L;
RH=5;%horizontal repitition if fixed

[~,dist2]=distcal1(spimg,targetimg,(bdctimg+Vglobal).*bdctsign,T);
fprintf('original distance %g\n\n',dist2);

for blockrow=1:RV
    for blockcol=1:RH
        %decide related coefficients
        if mod(NI,5)>=blockcol
            bwidth=floor(NI/5)+1;
        else
            bwidth=floor(NI/5);
        end        
        
        %calculate Cb
        tpm=tpm1(bdctimg+Vglobal,T,1);
        Cb=tpmtarget-tpm;
        Cb=Cb(:);
        fprintf('current Cb %g\n',norm(Cb));
        
        %initialize W        
        W=zeros(2*M+1,bwidth*L,(2*T+1)^2);
        vmask=true(2*M+1,bwidth*L);
        for row=1:L
            for col=1:bwidth
                if mod((blockrow-1)*L+row,8)==1 && mod((col-1)*5+blockcol,8)==1
                    vmask(:,(row-1)*bwidth+col)=false;
                    continue;
                end
                for f=1:2*M+1
                    flag=f-M-1;
                    if bdctimg((blockrow-1)*L+row,(col-1)*5+blockcol)+flag<0
                        vmask(f,(row-1)*bwidth+col)=false;
                        continue;
                    else
                        D=tpmdiff(bdctimg+Vglobal,(blockrow-1)*L+row,(col-1)*5+blockcol,flag,T);
                        W(f,(row-1)*bwidth+col,:)=reshape(D,1,1,(2*T+1)^2);
                    end
                end
            end
        end
        idxnode=find(vmask);
        lengthnode=length(idxnode);
        fprintf('size of network %d nodes\n',lengthnode);
        W=reshape(W,(2*M+1)*bwidth*L,(2*T+1)^2);
        W=W(vmask(:),:);
        %initalization of Tmat
        Tmat=zeros(lengthnode,lengthnode);
        
        %Since this matrix is symmetrical, we only need to calculate half of the
        %elements
        for m=1:lengthnode
            for n=m:lengthnode
                [x,i]=ind2sub([(2*M+1) (bwidth*L)],idxnode(m));
                [y,j]=ind2sub([(2*M+1) (bwidth*L)],idxnode(n));
                Tmat(m,n)=-A*delta(i,j)*(1-delta(x,y))-B-C*W(m,:)*W(n,:)';
            end
        end
        
        Tmat=Tmat+triu(Tmat,1)';
        
        %calculate I
        N=bwidth*L;
        I=N*B+C*W*Cb;
        
        %Calculate u00, U, V, and E
%         u00=u0*artanh(2*bwidth*L/lengthnode-1);
%         U=ones(lengthnode,1)*u00+(rand(lengthnode,1)*0.2-0.1)*u0;
        U=ones(2*M+1,bwidth*L)*(-1)*u0;
        U(M+1,:)=u0;
        U=U(vmask);
        V=nodeg(U,u0);
        E=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
        fprintf('iter:0 E0=%g\n',E);
        [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode,bwidth);
        fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n',fall,f1,f2,f3);
        
        %sequential update of nodes
        iter=0;
        while 1
            iter=iter+1;
            for m=1:lengthnode
                U(m)=Tmat(m,:)*V(:)+I(m);
                V(m)=nodeg(U(m),u0);
            end
            Enew=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
            fprintf('iter:%d  Enew=%g\n',iter,Enew);
            [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode,bwidth);
            fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n',fall,f1,f2,f3);
            stack=stafun((Enew-E)/abs(E),stack);
            %fprintf('length of stack %d\n',length(stack));
            %fprintf('mean of stack %g\n\n',mean(stack));
            if abs(mean(stack))<tol
                break
            end
            E=Enew;
        end
        
        %validation of V
        Vout=false(2*M+1,bwidth*L);
        Sout=repmat((-M:M)',1,bwidth*L);
        V=(V>0.5);
        Vout(idxnode)=V;
        Vout=sum(Vout.*Sout);
        
        newVglobal=Vglobal;
        newVglobal((blockrow-1)*L+1:(blockrow-1)*L+L,blockcol:5:(bwidth-1)*5+blockcol)=reshape(Vout,bwidth,L)';
        % result=networkvalidation(Vout);
        [~,dist2]=distcal1(spimg,targetimg,(bdctimg+newVglobal).*bdctsign,T);
        if dist2<norm(Cb)
            Vglobal=newVglobal;
            fprintf('current distance %g\n',dist2);
        else
            fprintf('new distance is worse than default, unchanged.\n');
        end
    end
end
bdctimg=(bdctimg+Vglobal).*bdctsign;

%objective function
function [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode,bwidth)
f1=0;
lengthnode=length(idxnode);
for m=1:lengthnode
    for n=1:lengthnode
        [x,i]=ind2sub([2*M+1 bwidth*L],idxnode(m));
        [y,j]=ind2sub([2*M+1 bwidth*L],idxnode(n));
        if i==j && x~=y
            f1=f1+V(m)*V(n);
        end
    end
end
f1=f1*A/2;

f2=B/2*(sum(V(:))-N)^2;

f3=repmat(V,[1 size(W,2)]).*W;
f3=sum(f3);
f3=sum((Cb-f3').^2);
f3=C/2*f3;
fall=f1+f2+f3;

%stack function: add newvalue into stack
function stack=stafun(newvalue,stack)
l=length(stack);
if l<3
    stack(l+1)=newvalue;
else
    stack(1:2)=stack(2:3);
    stack(3)=newvalue;
end

%ileagle output evaluation
% function result=networkvalidation(Vout)
% result=isequal(ones(1,size(Vout,2)),sum(Vout));
% if result==0
%     fprintf('not permutation matrix\n');
%     return
% end




    




    













