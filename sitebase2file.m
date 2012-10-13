function sitebase2file(fid,Dim,subnum,dim,gtypenum,gqn,ordermap)
%This function write the matlab generated site base information a to a
%binary file fid.

%Cheng Guo, 05 Aug 2009

fwrite(fid, Dim, 'int32');
fwrite(fid, subnum, 'int32');
fwrite(fid, dim, 'int32');
for i=1:subnum
    fwrite(fid, gtypenum, 'int32');
    for j=1:gtypenum
        fwrite(fid, gqn(j,i), 'int32');
    end
end
fwrite(fid, Dim, 'int32');
fwrite(fid, ordermap, 'int32');
