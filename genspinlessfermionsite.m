function site=genspinlessfermionsite(gqntype)
switch gqntype
    case 'fermion'
        %-spinless fermion base
        site.Dim=2;
        site.subnum=2;
        site.dim(1)=1;
        site.dim(2)=1;
        site.gtypenum=1;
        site.gqn(1,1)=0;
        site.gqn(1,2)=1;
        site.ordermap=zeros(site.Dim,1);
        
        
        %-Operators
        site.cm=[0 1;0 0];
        site.cp=[0 0;1 0];
        site.n=[0 0;0 1];
        site.sitedim=2;
    case 'spin'
end
end