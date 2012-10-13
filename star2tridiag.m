function [alpha,beta]=star2tridiag(indiag,inrow)

reps=0.01*min([min(indiag),min(inrow)]);
dim=length(indiag);

A=zeros(dim);

for n=1:dim
    A(n,n)=indiag(n);
    A(1,n)=inrow(n);
    A(n,1)=A(1,n);
end


alpha=zeros(dim,1);
beta=zeros(dim,1);

v=zeros(dim);
r=zeros(dim,1);
r(1)=1;
beta(1)=norm(r);
for j=1:dim
  % Basic recursion
  v(:,j)=r/beta(j);
  r=A*v(:,j);
  if j>1, r=r-v(:,j-1)*beta(j); end;
  alpha(j)=r'*v(:,j);
  r=r-v(:,j)*alpha(j);
  beta(j+1)=norm(r);
  % Reorthogonalize
  h1=v'*r;
  r=r-v*h1;
  nh1=norm(h1);
  if nh1>reps,
     r=r-v*(v'*r);
  end
end

alpha=alpha(2:end);
beta=beta(2:end-1);