function [bdctimg]=hopfieldnet(spimg,targetimg,T)
%Hopfieldnet: calculate the best BDCT with hopfield network
%intput images are square and multiple of 8

%Initialization
A=1500;
B=600;
C=40;
u0=300;
M=3;%maximum modification of coeff
flags=-M:-1;
flags=cat(2,flags,1:M);
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

for blockrow=1:RV
    for blockcol=1:RH
        %decide related coefficients
        if mod(NI,5)>=blockcol
            bwidth=floor(NI/5)+1;
        else
            bwidth=floor(NI/5);
        end        
        
        %calculate Cb
        %tpm=tpm1(abs((bdctimg+Vglobal).*bdctsign),T,1);
        tpm=tpm1(bdctimg+Vglobal,T,1);
        Cb=tpmtarget-tpm;
        Cb=Cb(:);
        fprintf('current Cb %g\n',norm(Cb));
        
        %initialize W        
        W=zeros(2*M+1,bwidth*L,(2*T+1)^2);
        vmask=true(2*M+1,bwidth*L);
        dcnum=0;
        for row=1:L
            for col=1:bwidth
                if mod((blockrow-1)*L+row,8)==1 && mod((col-1)*5+blockcol,8)==1
                    vmask(:,(row-1)*bwidth+col)=false;
                    dcnum=dcnum+1;
                    continue;
                end
                for f=flags
                    if bdctimg((blockrow-1)*L+row,(col-1)*5+blockcol)+f<0
                        vmask(f+M+1,(row-1)*bwidth+col)=false;
                        continue;
                    else
                        D=tpmdiff(bdctimg+Vglobal,(blockrow-1)*L+row,(col-1)*5+blockcol,f,T);
                        W(f+M+1,(row-1)*bwidth+col,:)=reshape(D,1,1,(2*T+1)^2);
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
        N=bwidth*L-dcnum;
        I=N*B+C*W*Cb;
        
        finish=false;
        tries=0;
        outoftry=false;
        while ~finish
            %Calculate u00, U, V, and E
            u00=u0*artanh(2*bwidth*L/lengthnode-1);
            U=ones(lengthnode,1)*u00+(rand(lengthnode,1)*0.2-0.1)*u0;
            %         U=ones(2*M+1,bwidth*L)*(-1)*u0;
            %         U(M+1,:)=u0;
            %         U=U(vmask);
            V=nodeg(U,u0);
            E=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
            fprintf('iter:0 E0=%g\n',E);
            %[fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode,bwidth);
            %fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n',fall,f1,f2,f3);
            
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
                %[fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode,bwidth);
                %fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n',fall,f1,f2,f3);
                stack=stafun(abs((Enew-E)/E),stack);
                %fprintf('length of stack %d\n',length(stack));
                %fprintf('mean of stack %g\n\n',mean(stack));
                if mean(stack)<tol
                    Vtest=V;
                    Vtest(Vtest>0.5)=1;
                    Vtest(Vtest<=0.5)=0;
                    [~,f1,f2,~]=objfun(A,B,C,Vtest,N,Cb,W,M,L,idxnode,bwidth);
                    if f1==0 && f2==0
                        break
                    end
                end
                E=Enew;
                if iter>30
                    outoftry=true;
                    break
                end
            end
            
            if outoftry
                continue
            end
            
            %validation of V
            Vout=zeros(2*M+1,bwidth*L);
            Sout=repmat((-M:M)',1,bwidth*L);
            V(V>0.5)=1;
            V(V<=0.5)=0;
            Vout(idxnode)=V;
            Vout=sum(Vout.*Sout);
            
            newVglobal=Vglobal;
            newVglobal((blockrow-1)*L+1:(blockrow-1)*L+L,blockcol:5:(bwidth-1)*5+blockcol)=reshape(Vout,bwidth,L)';
            newtpm=tpm1(bdctimg+newVglobal,T,1);
            dist2=norm(newtpm(:)-tpmtarget(:));
            if dist2<norm(Cb)
                Vglobal=newVglobal;
                fprintf('current distance %g\n',dist2);
                finish=true;
            else
                if tries<1
                    tries=tries+1;
                    fprintf('new distance is worse than default, reinitialising.\n');
                else
                    %use greedy algorithm to find solution
                    fprintf('using greedy algorithm\n');
                    newVglobal=Vglobal;
                    dist_old=norm(Cb);
                    
                    for row=1:L
                        for col=1:bwidth
                            if mod((blockrow-1)*L+row,8)==1 && mod((col-1)*5+blockcol,8)==1
                                continue
                            end
                            newtpm=tpm1(bdctimg+newVglobal,T,1);
                            flag=0;
                            for f=flags
                                if bdctimg((blockrow-1)*L+row,(col-1)*5+blockcol)+f<0
                                    continue
                                end
                                D=tpmdiff(bdctimg+newVglobal,(blockrow-1)*L+row,(col-1)*5+blockcol,f,T);
                                dist_new=norm(newtpm(:)+D(:)-tpmtarget(:));
                                if dist_new<dist_old
                                    flag=f;
                                    dist_old=dist_new;
                                end
                            end
                            if flag~=0
                                newVglobal((blockrow-1)*L+row,(col-1)*5+blockcol)=flag;
                            end
                        end
                    end
                    if dist_old>=norm(Cb)
                        fprintf('good solution not found by greedy algorithm.\n');
                        finish=true;
                    else
                        Vglobal=newVglobal;
                        fprintf('current distance is %g\n',dist_old);
                        finish=true;
                    end
                end
            end
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






    




    













