rdm=0;
tau=0.5;
T=100;
tarray=tau:tau:T;

folderlist=dir('Lambda*BathL*z*');
n=length(folderlist);
for k=1:n
filename=strcat(folderlist(k).name,'/results/rdm.dat');
if exist(filename)
    r=load(filename);
    rdm=rdm+r;
end
end

rdm=rdm/n;
plot(tarray,sqrt(rdm(:,1).^2+rdm(:,2).^2),'-o');