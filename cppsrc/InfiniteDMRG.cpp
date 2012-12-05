#include "InfiniteDMRG.h"

snake::physics::InfiniteDMRG::InfiniteDMRG()
{
}

void snake::physics::InfiniteDMRG::readParameters()
{
    std::cout<<"==========Reading model parameters from Matlab generated file."<<std::endl;
    std::string fname="./model/problemparmeters.dat";
    std::ifstream fin(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    fin.read((char*)&m_ChainLength,sizeof(int));
    int TGQN;
    fin.read((char*)&TGQN,sizeof(int));
    fin.close();
    TargetGQN[0]=TGQN;
    //TargetGQN.resize(3);
    //TargetGQN[1]=TargetGQN[0].gqn[0]+1;
    //TargetGQN[2]=TargetGQN[0].gqn[0]-1;

    fname="./model/Hfac.dat";
    fin.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    HoppingIntegrals.resize(m_ChainLength-1);
    OnSitePotentials.resize(m_ChainLength);
    TwoSitesInteraction.resize(m_ChainLength-1);
    snake::math::readvector(fin,HoppingIntegrals);
    snake::math::readvector(fin, OnSitePotentials);
    snake::math::readvector(fin, TwoSitesInteraction);
    fin.close();

    fname="./model/HC.dat";
    fin.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    snake::math::ReadOperators(fin,H,m_ChainLength-1);
    fin.close();

    std::cout<<"ChainLength="<<m_ChainLength<<std::endl;
    std::cout<<"TargetGQN="<<std::endl;
    snake::math::printvector(TargetGQN); std::cout<<std::endl;
    std::cout<<"Hopping Integrals are:"<<std::endl;
    snake::math::printvector(HoppingIntegrals);std::cout<<std::endl;
    std::cout<<"On site potentials are:"<<std::endl;
    snake::math::printvector(OnSitePotentials);std::cout<<std::endl;
    std::cout<<"Two sites nearest neighbor interaction are:"<<std::endl;
    snake::math::printvector(TwoSitesInteraction);std::cout<<std::endl;
}

void snake::physics::InfiniteDMRG::readFreesites()
{
    std::string fname="./model/site_operators.dat";
    std::ifstream opfin(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    std::string fname2="./model/site_base.dat";
    std::ifstream basefin(fname2.c_str(),std::ios_base::in|std::ios_base::binary);

    snake::physics::allfreesites.resize(m_ChainLength);
    for(int i=0;i<m_ChainLength;i++)
        snake::physics::allfreesites[i].readsite(basefin, opfin);
    opfin.close();
    basefin.close();
}


void snake::physics::InfiniteDMRG::run()
{
    readParameters();
    std::cout<<"==========Reading free site information from matlab generated files."<<std::endl;
    readFreesites();
    std::cout<<"==========Start Infinite DMRG."<<std::endl;

    ///The initial left and right block each with one site.
    left=new snake::physics::Chain(snake::physics::allfreesites[0],OnSitePotentials[0]);
    right=new snake::physics::Chain(snake::physics::allfreesites[m_ChainLength-1],OnSitePotentials[m_ChainLength-1]);
    left->write("./data/L");
    right->write("./data/R");
    if(m_ChainLength%2==1)
    {
        ///Add one more site to the initial left block if the ChainLengh is odd.
        ///This is used for the case when there is an extra impurity site at the left most.
        newleft=new snake::physics::Chain(*left,snake::physics::allfreesites[1],HoppingIntegrals[0], OnSitePotentials[1], TwoSitesInteraction[0]);
        delete left;
        left=newleft;
        left->write("./data/L");
    }
    int IniLeftL, IniRightL;
    IniLeftL=left->sitenum;
    IniRightL=right->sitenum;

    ///Determin the occupation number for the iDMRG
    std::vector<snake::math::GQN> tempTargetGQN;
    int tempTargetGQNNum=1;
    tempTargetGQN.resize(tempTargetGQNNum);
    double occnum;
    occnum=TargetGQN[0].gqn[0]/double(m_ChainLength);
    std::cout<<"The average occuptation number is "<<occnum<<std::endl;

    ///Increase both sides
    for(int i=1;i<=m_ChainLength/2-1;i++)
    {
        NewLeftL=IniLeftL+i;
        NewRightL=IniRightL+i;
            tempTargetGQN[0].gqn[0]=round((NewLeftL+NewRightL)*occnum);
            //std::cout<<tempTargetGQN[0].gqn[0]<<std::endl;
        if(i>1) delete m_SuperChain;
        addTwoSites(NewLeftL-1,NewRightL-1);
        m_SuperChain=new SuperChain(newleft,newright,left,right,H[NewLeftL-1]);
        m_SuperChain->TargetGQN=tempTargetGQN;
        //snake::math::printvector(tempTargetGQN);
        m_SuperChain->CalGroundState();
        //if(newleft->base.Dim>m_KeptStatesNum)
        m_SuperChain->renorm(m_KeptStatesNum);
        //else n++;
        newleft->write("./data/L");
        newright->write("./data/R");
        delete left;
        delete right;
        left=newleft;
        right=newright;
        std::cout<<std::endl;
    }
    delete newleft;
    delete newright;
    delete m_SuperChain;
    //std::cout<<std::endl<<"The maximium exact block size is "<<n<<std::endl<<std::endl;
    std::cout<<"==========Infinite DMRG complete."<<std::endl<<std::endl;
}


void snake::physics::InfiniteDMRG::addTwoSites(int LeftChainLength, int RightChainLength)
{
    newleft=new snake::physics::Chain(*left,snake::physics::allfreesites[LeftChainLength],HoppingIntegrals[LeftChainLength-1], OnSitePotentials[LeftChainLength], TwoSitesInteraction[LeftChainLength-1]);
    std::cout<<"NewLeftL="<<newleft->sitenum<<"\t";
    newright=new snake::physics::Chain(snake::physics::allfreesites[m_ChainLength-RightChainLength-1],
                                       *right,HoppingIntegrals[m_ChainLength-RightChainLength-1],
                                       OnSitePotentials[m_ChainLength-RightChainLength-1],
                                       TwoSitesInteraction[m_ChainLength-RightChainLength-1]);
    std::cout<<"NewRightL="<<newright->sitenum<<"\t";
}
