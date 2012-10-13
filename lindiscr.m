function [a,b] = lindiscr(Gamma, epsilon,HalfBandWidth)
% Input: the hybridization parameter Gamma and the energy level spacing epsilon
% Output: the Wilson chain hamiltonian diagnal term a and hopping term b.

% Based on Andreas's matlab program oliveira for logarithmic discretization

% Last modified 6 May 2008

  Edisc=HalfBandWidth:-epsilon:0;
  N=(length(Edisc)-1)*2;

% initialize Hamiltonian

  d=length(Edisc)-1; % number of intervalls
  c=sqrt(2*Gamma/pi);

  d0=0.5*(Edisc(1:end-1)+ Edisc(2:end)); % diagonal entries
  d1=c/sqrt(2) * sqrt(abs(diff(Edisc))); % coupling to impurity

  hdd=sparse(diag(d0)); Z=sparse(d,d);

  H0=[ 0     d1    d1
       d1'  +hdd   Z
       d1'   Z    -hdd ];

  D=size(H0,1);    % 2*ndim + 1
  niter=D;         % max number of iterations


% starting lanczos tridiagonalization
  U=zeros(D,1); U(1)=1; % start vector

  alpha=zeros(1,niter); beta=zeros(1,niter); betaMin=1E-16;

  for i=1:niter
     v=H0*U(:,i);
     alpha(i)=U(:,i)'*v;

     v=v-U*(U'*v); v=v-U*(U'*v); % twice for numerical reasons
     beta(i)=norm(v);

     if (beta(i)<betaMin)
        break
     else U(:,i+1)=v/beta(i); end
  end

  alpha=alpha(1:i); beta=beta(1:i-1);
  a=alpha; b=beta;

return


