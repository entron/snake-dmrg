
tau=0.5;
T=2000;
tarray=tau:tau:T;
numofz=4;

folderlist=dir('ratio*Lambda1.1BathL*');
n=length(folderlist);
sigmazfinal=zeros(n/numofz,1);
omega=sigmazfinal;
index=1;
for k=1:n
    k
    sigmaz=0;
    m=k;
 %   for m=k:k+numofz-1
        filename=strcat(folderlist(m).name,'/results/sigmaz_t.dat');
        if exist(filename)
            r=load(filename);
            sigmaz=sigmaz+r;
            %sscanf(filename,'ratio%g*')
        end
  %  end
    omega(index)=sscanf(filename,'ratio%g*');
    sigmazlast=sigmaz(end-100:end);
    [maxval,maxindex]=max(sigmazlast);
    [minval,minindex]=min(sigmazlast);
    if maxindex>minindex
        a=minindex;b=maxindex;
    else
        a=maxindex;b=minindex;
    end
    %sigmazfinal(index)=mean(sigmazlast(a:b));
    sigmazfinal(index)=(maxval+minval)/2;
    index=index+1;
end


p=(1+sigmazfinal)./2;
p=cat(2,omega,p);
p=sortrows(p);
plot(p(:,1),p(:,2),'-o');