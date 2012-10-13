function savepara(para)

%%%%%%%%%%%%%%%%%%-------WRITE TO FILES--------%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%------Problem Parameters
filename=[para.folder '/problemparmeters.dat'];
fid=fopen(filename, 'wb');
fwrite(fid,para.L,'int32');
fwrite(fid,para.TGQN,'int32');
fwrite(fid,para.t_num,'int32');
fclose(fid);

%%%%%%------Site info
filename=[para.folder '/site_base.dat'];
fid=fopen(filename, 'wb');
sitebase2file(fid,para.site{1}.Dim,para.site{1}.subnum,para.site{1}.dim,para.site{1}.gtypenum,para.site{1}.gqn,para.site{1}.ordermap);
for i=2:para.L
sitebase2file(fid,para.site{2}.Dim,para.site{2}.subnum,para.site{2}.dim,para.site{2}.gtypenum,para.site{2}.gqn,para.site{2}.ordermap);
end
fclose(fid);

filename=[para.folder '/site_operators.dat'];
fid=fopen(filename, 'wb');
mat2file(fid,para.site{1}.sigmax);
mat2file(fid,para.site{1}.sigmaz);
for i=2:para.L
    mat2file(fid,para.site{2}.cm);
    mat2file(fid,para.site{2}.n);
end
fclose(fid);

%%%%%%-----Starting State Hamiltonian
filename=[para.folder '/Hfac.dat'];
fid=fopen(filename, 'wb');
fwrite(fid,para.hopT, 'real*8');
fwrite(fid,para.onesiteE,'real*8');
fwrite(fid,para.twositesV,'real*8');
fclose(fid);

filename=[para.folder '/HC.dat'];
fid=fopen(filename, 'wb');
for i=1:para.L-1
    mat2file(fid, para.h0{i});
end
fclose(fid);

%%%%%%-----Time evolving Hamiltonian
filename=[para.folder '/rt_T0.dat'];
fid=fopen(filename, 'wb');
for i=1:para.L-1
mat2file(fid,para.U{i});
end
fclose(fid);

filename=[para.folder '/rt_H1_T0.dat'];
fid=fopen(filename, 'wb');
for i=1:para.t_num
	mat2file(fid, para.U_imp{i});
end
fclose(fid);
end
