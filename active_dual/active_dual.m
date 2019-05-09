function [x_opt,u_opt]=active_dual(Q,R,xmin,xmax,umin,umax,T,x0,A,B,w,xf)
n = size(Q,1);
%m = size(R,1);
invQ=inv(2*Q);
invR=inv(2*R);
la0=ones((T+1)*n,1);
la=la0;
[x,u,indx,indu]=find_zstar(la0,invQ,invR,A,B,xmin,xmax,umin,umax,T);

Dla=calderi(x,u,A,B,w,x0,xf,T);
Hla=calhess(indx,indu,invQ,invR,A,B,T);
dla=linsolve(Hla,-Dla);
t=1;
la1=la+t*dla;
f1=calf(la1,Q,R,A,B,xmin,xmax,umin,umax,T,[],[],x0,xf,w);
f0=calf(la,Q,R,A,B,xmin,xmax,umin,umax,T,x,u,x0,xf,w);
while f1<f0+0.01*t*Dla'*dla
    t=0.5*t;
    la1=la+t*dla;
    f1=calf(la1,Q,R,A,B,xmin,xmax,umin,umax,T,[],[],x0,xf,w);
end
la=la1;
[x,u,indx,indu]=find_zstar(la,invQ,invR,A,B,xmin,xmax,umin,umax,T);
Dla=calderi(x,u,A,B,w,x0,xf,T);
num=1;
while norm(Dla)>1e-2
Hla=calhess(indx,indu,invQ,invR,A,B,T);

dla=linsolve(Hla,-Dla);
t=1;

la1=la+t*dla;

f1=calf(la1,Q,R,A,B,xmin,xmax,umin,umax,T,[],[],x0,xf,w);
f0=calf(la,Q,R,A,B,xmin,xmax,umin,umax,T,x,u,x0,xf,w);

 while f1<f0+0.01*t*Dla'*dla
     t=0.5*t;
     if t<0.00001
         break
     end
     la1=la+t*dla;
     f1=calf(la1,Q,R,A,B,xmin,xmax,umin,umax,T,[],[],x0,xf,w);
 end

la=la1;

[x,u,indx,indu]=find_zstar(la,invQ,invR,A,B,xmin,xmax,umin,umax,T);
Dla=calderi(x,u,A,B,w,x0,xf,T);
num=num+1;
 if num>1e4
     break
 end


end
x_opt=x;
u_opt=u;
end


function [x0,u0,xindd,uindd]=find_zstar(la,invQ,invR,A,B,xmin,xmax,umin,umax,T)
n=size(invQ,1);
m=size(invR,1);
u_all=zeros(m*T,3);
x_all=zeros(n*T,3);

u_all(:,1)=repmat(umax,T,1);
u_all(:,3)=repmat(umin,T,1);
x_all(:,1)=repmat(xmax,T,1);
x_all(:,3)=repmat(xmin,T,1);
bdB=kron(eye(T),B');
bdinvR=kron(eye(T),invR);
u_all(:,2)=-bdinvR*bdB*la(1:n*T);
x_all(n*(T-1)+1:n*T,2)=-invQ*(la(T*n+1:(T+1)*n)-la((T-1)*n+1:T*n));
bdA=kron(eye(T-1),A');
bdinvQ=kron(eye(T-1),invQ);
x_all(1:n*(T-1),2)=-bdinvQ*(bdA*la(n+1:T*n)-la(1:n*(T-1)));
[~,xindd]=max(x_all,[],2);
[~,uindd]=max(u_all,[],2);

for i=1:n*T
x_all(i,xindd(i))=-inf;
end
for i=1:m*T
u_all(i,uindd(i))=-inf;
end
[x0,xindd]=max(x_all,[],2);

[u0,uindd]=max(u_all,[],2);

end

function Dla=calderi(x,u,A,B,w,x0,xf,T)
n=size(A,2);
m=size(B,2);
D_dis=zeros((T+1)*n,T+1);
D_dis(1:n,1)=A*x0+B*u(1:m)+w;
for i=2:T
    D_dis((i-2)*n+1:(i-1)*n,i)=-x((i-2)*n+1:(i-1)*n);
    D_dis((i-1)*n+1:i*n,i)=A*x((i-2)*n+1:(i-1)*n)+B*u((i-1)*m+1:i*m)+w;
end
D_dis((T-1)*n+1:T*n,T+1)=-x((T-1)*n+1:T*n);
D_dis(T*n+1:(T+1)*n,T+1)=x((T-1)*n+1:T*n)-xf;
Dla=sum(D_dis,2);
end

function Hla=calhess(indx,indu,invQ,invR,A,B,T)
n=size(invQ,1);
m=size(invR,1);
Hla=zeros((T+1)*n,(T+1)*n);
x_ind=zeros(n,T);
u_ind=zeros(m,T);
AA=zeros(n,n,n);
BB=zeros(n,n,m);

for i=1:n
    AA(:,:,i)=A(:,i)*(A(:,i))';
end

for i=1:m
    BB(:,:,i)=B(:,i)*(B(:,i))';
end

for i=1:T
    indd=find(indx((i-1)*n+1:i*n)==2);
    x_ind(indd,i)=1;
    indd=find(indu((i-1)*m+1:i*m)==2);
    u_ind(indd,i)=1;
end

Hla(1:n,1:n)=-diag(x_ind(:,1))*invQ;


for i=1:m
    Hla(1:n,1:n)=Hla(1:n,1:n)-u_ind(i,1)*invR(i,i)*BB(:,:,i);
end

Hla(1:n,n+1:2*n)=diag(x_ind(:,1))*invQ*A;
for i=2:T
    Hla((i-1)*n+1:i*n,(i-2)*n+1:(i-1)*n)= Hla((i-2)*n+1:(i-1)*n,(i-1)*n+1:i*n)';
    Hla((i-1)*n+1:i*n,(i-1)*n+1:i*n)=-diag(x_ind(:,i))*invQ;
    for j=1:m
    Hla((i-1)*n+1:i*n,(i-1)*n+1:i*n)=Hla((i-1)*n+1:i*n,(i-1)*n+1:i*n)-u_ind(j,i)*invR(j,j)*BB(:,:,j);
    end
    for j=1:n
    Hla((i-1)*n+1:i*n,(i-1)*n+1:i*n)=Hla((i-1)*n+1:i*n,(i-1)*n+1:i*n)-x_ind(j,i-1)*invQ(j,j)*AA(:,:,j);
    end
    
    Hla((i-1)*n+1:i*n,i*n+1:(i+1)*n)=diag(x_ind(:,i))*invQ*A;
end

Hla(T*n+1:(T+1)*n,T*n+1:(T+1)*n)=-diag(x_ind(:,T))*invQ;


Hla(T*n+1:(T+1)*n,(T-1)*n+1:T*n)=Hla((T-1)*n+1:T*n,T*n+1:(T+1)*n)';
end

function  f=calf(la,Q,R,A,B,xmin,xmax,umin,umax,T,x,u,x0,xf,w)
n=size(Q,1);
m=size(R,1);
if isempty(x)&& isempty(u)
    invQ=inv(2*Q);
    invR=inv(2*R);
    [x,u,~,~]=find_zstar(la,invQ,invR,A,B,xmin,xmax,umin,umax,T);
end

 bdQ=kron(eye(T),Q);
 bdR=kron(eye(T),R);
 bdA=kron(eye(T-1),A);
 bdB=kron(eye(T-1),B);
 f=x'*bdQ*x+u'*bdR*u+la(1:n)'*(A*x0+B*u(1:m)+w-x(1:n))+la(n+1:T*n)'*(bdA*x(1:(T-1)*n)+bdB*u(m+1:T*m)-x(n+1:T*n)+repmat(w,T-1,1))+la(n*T+1:(T+1)*n)'*(x((T-1)*n+1:T*n)-xf);
end