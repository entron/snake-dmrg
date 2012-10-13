function mat2file(fid,a)
%This function write a matlab matrix a to a binary file fid according to
%fortran storage type.

%Cheng Guo, 15 Dec 2008

[row,col]=size(a);
fwrite(fid,row,'int32');
fwrite(fid,col,'int32');

if isreal(a)
    fwrite(fid,a,'real*8');
else
    for i=1:col
        for j=1:row
            fwrite(fid,real(a(j,i)),'real*8');
            fwrite(fid,imag(a(j,i)),'real*8');
        end
    end
end
