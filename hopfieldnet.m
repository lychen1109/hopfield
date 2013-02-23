function [Tmat,I]=hopfieldnet(spimg,targetimg)
%Hopfieldnet: calculate the best BDCT with hopfield network
%expected output bdctimg
%todos: 1) filter out those illegal Vxi; 
%       2) add input params which can shift colm;
%       3) test if T can be more than 4
%       4) add input for M according to memory of machine
%       5) compare different choices of A, B, C
%       6) test single precision storage for Tmat

%Initialization
A=1500;
B=400;
C=2;
u0=300;
T=4;%Threshold of markov algorithm
M=3;%maximum modification of coeff
L=3;%Process L lines a time
N=25*L+5;
stack=[];%used in calculation of average delta Energy function
tol=5e-4;%tolerance of minimum delta Energy function
fprintf('size of network %d nodes\n',(2*M+1)*25*L);

%calculate tpm of spimg and targetimg
spimg=reshape(spimg,128,128);
targetimg=reshape(targetimg,128,128);
bdctimg=blkproc(spimg,[8 8],@dct2);
bdctsign=sign(bdctimg);
bdctimg=abs(round(bdctimg));
bdcttarget=blkproc(targetimg,[8 8],@dct2);
bdcttarget=abs(round(bdcttarget));
tpm=tpm1(bdctimg,T,1);
tpmtarget=tpm1(bdcttarget,T,1);

%calculate Cb
Cb=tpmtarget-tpm;
Cb=Cb(:);

%initialize W
W=zeros(2*M+1,25*L,(2*T+1)^2);
for row=1:L
    for col=1:25
        for f=1:2*M+1
            flag=f-M-1;
            if bdctimg(row,(col-1)*5+3)+flag<0
                continue;
            else
                D=tpmdiff(bdctimg,row,(col-1)*5+3,flag);
            end
            W(f,(row-1)*25+col,:)=reshape(D,1,1,(2*T+1)^2);
        end
    end
end

%initalization of Tmat
Tmat=zeros((2*M+1)*25*L,(2*M+1)*25*L);

%Since this matrix is symmetrical, we only need to calculate half of the
%elements
for m=1:(2*M+1)*25*L
    for n=m:(2*M+1)*25*L
        [x,i]=ind2sub([(2*M+1) (25*L)],m);
        [y,j]=ind2sub([(2*M+1) (25*L)],n);
        Tmat(m,n)=-A*delta(i,j)*(1-delta(x,y))-B-C*dot(W(x,i,:),W(y,j,:),3);
    end
end

Tmat=Tmat+triu(Tmat,1)';

%calculate I
I=ones(2*M+1,25*L)*N*B;
for x=1:2*M+1
    for i=1:25*L
        I(x,i)=I(x,i)+C*(Cb'*squeeze(W(x,i,:)));
    end
end

%Calculate u00, U, V, and E
u00=u0*artanh(2/(2*M+1)-1);
U=ones(2*M+1,25*L)*u00+(rand(2*M+1,25*L)*0.2-0.1)*u0;
V=nodeg(U,u0);
E=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
fprintf('E0=%g\n',E);
[fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W);
fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n\n',fall,f1,f2,f3);

%sequential update of nodes
iter=0;
while 1
    iter=iter+1;
    for x=1:2*M+1
        for i=1:25*L
            m=sub2ind([2*M+1 25*L],x,i);
            U(x,i)=Tmat(m,:)*V(:)+I(m);
            V(x,i)=nodeg(U(x,i),u0);
        end
    end
    Enew=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
    fprintf('iter:%d  Enew=%g\n',iter,Enew);
    [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W);
    fprintf('fall=%g, f1=%g, f2=%g, f3=%g\n',fall,f1,f2,f3);
    stack=stafun((Enew-E)/abs(E),stack);
    %fprintf('length of stack %d\n',length(stack));
    %fprintf('mean of stack %g\n\n',mean(stack));
    if abs(mean(stack))<tol
        break
    end
    E=Enew;
end

%objective function
function [fall,f1,f2,f3]=objfun(A,B,C,V,N,Cb,W)
[sx,si]=size(V);
f1=0;
for i=1:si
    for x=1:sx-1
        for y=x+1:sx
            f1=f1+V(x,i)*V(y,i);
        end
    end
end
f1=f1*A;%No need to be divided by 2

f2=B/2*(sum(V(:))-N)^2;

f3=repmat(V,[1 1 size(W,3)]).*W;
f3=sum(f3,1);
f3=sum(f3,2);
f3=sum((Cb-squeeze(f3)).^2);
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




    













