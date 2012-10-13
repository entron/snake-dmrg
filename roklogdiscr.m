function [alpha,beta,N]=roklogdiscr(Gamma,D,Lambda,BathL,z)

%Calculate the wilson chain parameters for spinless Fermion model
%See Appendix A of PRB 71, 045122 for details.
%Gamma is the hybridizaiton parameter, D is the half band width, Lambda is the logarithmic discretization parameter and N is the Wilson chain length.

J=Gamma;
tempfac=(1-Lambda^(-1))/log(Lambda);
N=BathL/2;
xi_positive=zeros(N,1);
gamma_positive=zeros(N,1);

w0=D;
wdata=zeros(N+1,1);
wdata(1)=w0;
for j=2:N+1
    w1=D*Lambda.^(2-j-z);
    wdata(j)=w1;
    gamma_positive(j-1)=sqrt((w0-w1)*J);
    if j==2
	    xi_positive(j-1)=D*((1-Lambda^(-1))/log(Lambda)+1-z);
    else
    	    xi_positive(j-1)=D*tempfac*Lambda^(3-j-z);
    end
    w0=w1;
end

fprintf('wdata(end-5:end) = %.10g\n', wdata(end-5:end));
%bar(log(wdata),ones(N+1,1));
% gamma(end-5:end)
% xi(end-5:end)



indiag=zeros(BathL+1,1);
inrow=indiag;

for n=1:N
    inrow(n+1)=gamma_positive(n)/sqrt(pi);
    inrow(BathL-n+2)=inrow(n+1);
    indiag(n+1)=xi_positive(n);
    indiag(BathL-n+2)=-1*indiag(n+1);
end
% indiag(end-5:end)
% inrow(end-5:end)
[alpha,beta]=star2tridiag(indiag,inrow);
