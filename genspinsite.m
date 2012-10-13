function site=genspinsite(gqntype)
switch gqntype
    case 'fermion' %Therefore subnum=1 and do not use spin as good quantum number.
        %-spin base
        site.Dim=2;
        site.subnum=1;
        site.dim=2;
        site.gtypenum=1;
        site.gqn(1,1)=0;
        site.ordermap=zeros(site.Dim,1);
        
        %-spin operators
        site.sigmax=[0 1;1 0];
        site.sigmay=[0 -1i;1i 0];
        site.sigmaz=[-1 0;0 1];
        site.sitedim=2;
    case 'spin'
end
end