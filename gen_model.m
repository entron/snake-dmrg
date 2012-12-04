function gen_model(Omega,omegaratio,Lambda,BathL,z)


%%%%%%%%%%%%%%%%%%---------MODEL INFO-------%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%-------Model Parameters
FlorianGamma=0.1;
D=1;
Delta=0.2;
vz=0;
vx=0.1;
epsilon0=0;
ed0=100;
%Driven parameters
w=omegaratio*Delta; %driving frequesncy

%%%%%------Discretization Parameters
%Lambda=2;
%BathL=20;
inter=D*2/BathL; [a,b]=lindiscr(FlorianGamma/2,inter,D);
%[a,b]=roklogdiscr(FlorianGamma/2,D,Lambda,BathL,z);
para.TGQN=BathL/2
para.L=BathL+2

folder=sprintf('ratio%gLambda%gBathL%gz%g',omegaratio,Lambda,BathL,z);
mkdir(folder);
system(['cp ./build/cppsrc/Snake ' folder]);
para.folder=[folder,'/model'];
mkdir(para.folder);


%%%%%------Time evolution Parameters
%%Trotter step of real time
rt_tau=0.5;
%Stating time
t_start=0;
%Ending time
t_end=2000;
%Time span
T=t_end-t_start;
%Time step number
para.t_num=T/rt_tau;


%%%%%%%%%%%%%%%%%%---------SITE INFO--------%%%%%%%%%%%%%%%%%%%%%%%%%
para.site{1}=genspinsite('fermion');
para.site{2}=genspinlessfermionsite('fermion');


%%%%%%%%%%%%%%------STARTING STATE HAMILTONIAN-------%%%%%%%%%%%%%%%%%
para.onesiteE=zeros(para.L,1);
para.onesiteE(1)=ed0;
para.twositesV=zeros(para.L-1,1);
para.hopT=zeros(para.L-1,1);
para.hopT(1)=0;
for i=2:para.L-1
    para.hopT(i)=b(i-1);
end

%para.h0{1}=ed0*kron(para.site{1}.sigmax,eye(2));
para.h0first=ed0*para.site{1}.sigmaz;
para.h0{1}=kron(para.h0first,eye(2));
for i=2:para.L-1
para.h0{i}=para.hopT(i)*(kron(para.site{2}.cm,para.site{2}.cp)+kron(para.site{2}.cp,para.site{2}.cm));
end 


%%%%%%%%%%%%%%------TIME EVOLVING HAMILTONIAN--------%%%%%%%%%%%%%%%%%

%%%%%%-------constat Wilson chain part 
para.h{1}=eye(para.site{1}.sitedim*para.site{2}.sitedim); %Useless term
para.h{2}=epsilon0*kron(para.site{2}.n,eye(2))+para.hopT(2)*(kron(para.site{2}.cm,para.site{2}.cp)+kron(para.site{2}.cp,para.site{2}.cm));
for i=3:para.L-1
para.h{i}=para.hopT(i)*(kron(para.site{2}.cm,para.site{2}.cp)+kron(para.site{2}.cp,para.site{2}.cm));
end 

for i=1:para.L-1
para.U{i}=expm(-1i*(rt_tau/2)*para.h{i});
end

%%%%%%-------Time-dependent impurity term
for i=1:para.t_num
    t=(i-0.5)*rt_tau+t_start;
    %para.h_imp_t{i}=(Delta/2)*kron(para.site{1}.sigmaz,eye(2))+(vz/2)*kron(para.site{1}.sigmaz,para.site{2}.n)+(vx/2)*kron(para.site{1}.sigmax,para.site{2}.n);
    para.h_imp_t{i}=Omega*cos(w*t)*kron(para.site{1}.sigmax,eye(2))+(Delta/2)*kron(para.site{1}.sigmaz,eye(2))+(vz/2)*kron(para.site{1}.sigmaz,para.site{2}.n)+(vx/2)*kron(para.site{1}.sigmax,para.site{2}.n);
    para.U_imp{i}=expm(-1i*(rt_tau/2)*para.h_imp_t{i});
end

%Save parameters to file for c++ program to read
savepara(para)

end






