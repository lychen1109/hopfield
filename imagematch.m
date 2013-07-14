function V=imagematch(dismat,opt)
%match image based on their distance matric
% the output of the function is a connection matrix of two image list

%initialization
A=opt.A;
B=opt.B;
C=opt.C;
D=opt.D;
u0=opt.u0;
Tol=1e-5;

%calculate Tmat and I
T11=ones(90,90);
T11=T11-eye(90);
diagcell=cell(1);
diagcell{1}=T11;
diagcell=repmat(diagcell,90,1);
T1=blkdiag(diagcell{:});
T2=zeros(90^2,90^2);
for i=1:89
    T2=T2+diag(ones(90^2-i*90,1),i*90);
end
T2=T2+T2';
T3=ones(90^2,90^2);
T4=dismat(:)*dismat(:)';
Tmat=-A*T1-B*T2-C*T3-D*T4;
N=90;
I=N*C*ones(90^2,1);

%initialize U
u00=u0*artanh(-44/45);
U=ones(90^2,1)*u00+(rand(90^2,1)*0.2-0.1)*u0;
V=nodeg(U,u0);
E=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
fprintf('iter:0  E0=%g\n',E);
[f1,f2,f3,f4]=objfun(A,B,C,D,V,N,dismat,T1,T2);
fprintf('f1=%g, f2=%g, f3=%g, f4=%g, f=%g\n',f1,f2,f3,f4,f1+f2+f3+f4);

iter=0;
while 1
    iter=iter+1;
    idx=randperm(90^2);
    for i=1:90^2
        m=idx(i);
        U(m)=Tmat(m,:)*V(:)+I(m);
        V(m)=nodeg(U(m),u0);
    end
    Enew=-0.5*V(:)'*Tmat*V(:)-V(:)'*I(:);
    fprintf('iter:%d  Enew=%g\n',iter,Enew);
    [f1,f2,f3,f4]=objfun(A,B,C,D,V,N,dismat,T1,T2);
    fprintf('f1=%g, f2=%g, f3=%g, f4=%g, f=%g\n',f1,f2,f3,f4,f1+f2+f3+f4);
    delta=(Enew-E)/E;
    fprintf('delta E=%g\n\n',delta);
    if abs(delta)<Tol
        break;
    else
        E=Enew;
    end
end

V=reshape(V,90,90);
V(V>0.5)=1;
V(V<=0.5)=0;
if ~isequal(sum(V),ones(1,90))
    fprintf('Not permutation matrix\n');
end

function [f1,f2,f3,f4]=objfun(A,B,C,D,V,N,dismat,T1,T2)
crossmat=V*V';
f1=A/2*sum(sum(crossmat.*T1));
f2=B/2*sum(sum(crossmat.*T2));
f3=C/2*(sum(V(:))-N)^2;
f4=D/2*sum(V(:).*dismat(:))^2;







