function V=hopfieldnet(spimg,targetimg)
%Hopfieldnet: calculate the best BDCT with hopfield network
%expected output bdctimg
%todos: 1) filter out those illegal Vxi; 
%       2) add input params which can shift colm;
%       3) test if T can be more than 4
%       4) add input for M according to memory of machine
%       5) compare different choices of A, B, C
%       6) test single precision storage for Tmat

%Initialization
A=500;
B=200;
C=500;
u0=0.02;
T=4;%Threshold of markov algorithm
M=1;%maximum modification of coeff
L=1;%Process L lines a time
N=25*L+4;

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
U=ones(2*M+1,25*L)*u00;
V=nodeg(U,u0);
E=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
fprintf('E0=%d\n',E);

%sequential update of nodes
for loop=1:20
    for x=1:2*M+1
        for i=1:25*L
            m=sub2ind([2*M+1 25*L],x,i);
            U(x,i)=Tmat(m,:)*V(:)+I(m);
            V(x,i)=nodeg(U(x,i),u0);
        end
    end    
    Enew=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
    fprintf('Enew=%d\n',Enew);
    if Enew>=E
        break
    else
        E=Enew;
    end
end













