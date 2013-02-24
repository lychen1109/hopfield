function [Vout]=hopfieldnet(bdctimg,tpmtarget,T)
%Hopfieldnet: calculate the best BDCT with hopfield network
%expected output bdctimg
%todos: 
%       2) add input params which can shift colm;

%Initialization
A=1500;
B=400;
C=2;
u0=300;
M=3;%maximum modification of coeff
L=3;%Process L lines a time
N=25*L+5;
stack=[];%used in calculation of average delta Energy function
tol=5e-4;%tolerance of minimum delta Energy function

%calculate current tpm of spimg
tpm=tpm1(bdctimg,T,1);

%calculate Cb
Cb=tpmtarget-tpm;
Cb=Cb(:);

%initialize W
W=zeros(2*M+1,25*L,(2*T+1)^2);
vmask=true(2*M+1,25*L);
for row=1:L
    for col=1:25
        for f=1:2*M+1
            flag=f-M-1;
            if bdctimg(row,(col-1)*5+3)+flag<0
                vmask(f,(row-1)*25+col)=false;
                continue;
            else
                D=tpmdiff(bdctimg,row,(col-1)*5+3,flag,T);
                W(f,(row-1)*25+col,:)=reshape(D,1,1,(2*T+1)^2);
            end
        end
    end
end
idxnode=find(vmask);
lengthnode=length(idxnode);
fprintf('size of network %d nodes\n',lengthnode);
W=reshape(W,(2*M+1)*25*L,(2*T+1)^2);
W=W(vmask(:),:);
%initalization of Tmat
Tmat=zeros(lengthnode,lengthnode);

%Since this matrix is symmetrical, we only need to calculate half of the
%elements
for m=1:lengthnode
    for n=m:lengthnode
        [x,i]=ind2sub([(2*M+1) (25*L)],idxnode(m));
        [y,j]=ind2sub([(2*M+1) (25*L)],idxnode(n));
        Tmat(m,n)=-A*delta(i,j)*(1-delta(x,y))-B-C*W(m,:)*W(n,:)';
    end
end

Tmat=Tmat+triu(Tmat,1)';

%calculate I
I=N*B+C*W*Cb;

%Calculate u00, U, V, and E
u00=u0*artanh(2*25*L/lengthnode-1);
U=ones(lengthnode,1)*u00+(rand(lengthnode,1)*0.2-0.1)*u0;
V=nodeg(U,u0);
E=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
fprintf('E0=%g\n',E);
% [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode);
% fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n\n',fall,f1,f2,f3);

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
%     [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode);
%     fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n',fall,f1,f2,f3);
    stack=stafun((Enew-E)/abs(E),stack);
    %fprintf('length of stack %d\n',length(stack));
    %fprintf('mean of stack %g\n\n',mean(stack));
    if abs(mean(stack))<tol
        break
    end
    E=Enew;
end

%validatin of V
Vout=false(2*M+1,25*L);
V=(V>0.5);
Vout(idxnode)=V;
% result=networkvalidation(Vout);


%objective function
% function [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W,M,L,idxnode)
% f1=0;
% lengthnode=length(idxnode);
% for m=1:lengthnode
%     for n=1:lengthnode
%         [x,i]=ind2sub([2*M+1 25*L],idxnode(m));
%         [y,j]=ind2sub([2*M+1 25*L],idxnode(n));
%         if i==j && x~=y
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




    




    













