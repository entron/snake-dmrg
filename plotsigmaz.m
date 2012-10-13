sigmaz=0;
tau=0.5;
T=2000;
tarray=tau:tau:T;

folderlist=dir('ratio1.68*');
n=length(folderlist);
for k=1:n
filename=strcat(folderlist(k).name,'/results/sigmaz_t.dat');
if exist(filename)
    r=load(filename);
    sigmaz=sigmaz+r;
end
end

sigmaz=sigmaz/n;
plot(tarray,sigmaz,'-o');