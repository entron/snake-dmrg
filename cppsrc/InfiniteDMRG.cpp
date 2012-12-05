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
    m_TargetGQN[0]=TGQN;

    fname="./model/Hfac.dat";
    fin.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    m_HoppingIntegrals.resize(m_ChainLength-1);
    m_OnSitePotentials.resize(m_ChainLength);
    m_TwoSitesInteraction.resize(m_ChainLength-1);
    snake::math::readvector(fin,m_HoppingIntegrals);
    snake::math::readvector(fin, m_OnSitePotentials);
    snake::math::readvector(fin, m_TwoSitesInteraction);
    fin.close();

    fname="./model/HC.dat";
    fin.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    snake::math::ReadOperators(fin,m_TwoFreeSitesHamiltonian,m_ChainLength-1);
    fin.close();

    printParameters();
}

void snake::physics::InfiniteDMRG::printParameters()
{
    std::cout<<"ChainLength="<<m_ChainLength<<std::endl;
    std::cout<<"TargetGQN="<<std::endl;
    snake::math::printvector(m_TargetGQN); std::cout<<std::endl;
    std::cout<<"Hopping Integrals are:"<<std::endl;
    snake::math::printvector(m_HoppingIntegrals);std::cout<<std::endl;
    std::cout<<"On site potentials are:"<<std::endl;
    snake::math::printvector(m_OnSitePotentials);std::cout<<std::endl;
    std::cout<<"Two sites nearest neighbor interaction are:"<<std::endl;
    snake::math::printvector(m_TwoSitesInteraction);std::cout<<std::endl;
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

void snake::physics::InfiniteDMRG::initLeftRightChains()
{
    ///The initial left and right block each with one site.
    m_LeftChain=new snake::physics::Chain(snake::physics::allfreesites[0],m_OnSitePotentials[0]);
    m_RightChain=new snake::physics::Chain(snake::physics::allfreesites[m_ChainLength-1],m_OnSitePotentials[m_ChainLength-1]);
    m_LeftChain->write("./data/L");
    m_RightChain->write("./data/R");
    if(m_ChainLength%2==1)
    {
        ///Add one more site to the initial left block if the ChainLengh is odd.
        ///This is used for the case when there is an extra impurity site at the left most.
        m_NewLeftChain=new snake::physics::Chain(*m_LeftChain,snake::physics::allfreesites[1],m_HoppingIntegrals[0], m_OnSitePotentials[1], m_TwoSitesInteraction[0]);
        delete m_LeftChain;
        m_LeftChain=m_NewLeftChain;
        m_LeftChain->write("./data/L");
    }
}


void snake::physics::InfiniteDMRG::run()
{
    readParameters();
    std::cout<<"==========Reading free site information from matlab generated files."<<std::endl;
    readFreesites();
    std::cout<<"==========Start Infinite DMRG."<<std::endl;

    initLeftRightChains();

    int IniLeftL, IniRightL;
    IniLeftL=m_LeftChain->sitenum;
    IniRightL=m_RightChain->sitenum;

    ///Determin the occupation number for the iDMRG
    std::vector<snake::math::GQN> tempTargetGQN;
    int tempTargetGQNNum=1;
    tempTargetGQN.resize(tempTargetGQNNum);
    double occnum;
    occnum=m_TargetGQN[0].gqn[0]/double(m_ChainLength);
    std::cout<<"The average occuptation number is "<<occnum<<std::endl;

    ///Increase both sides
    for(int i=1;i<=m_ChainLength/2-1;i++)
    {
        m_NewLeftL=IniLeftL+i;
        m_NewRightL=IniRightL+i;
            tempTargetGQN[0].gqn[0]=round((m_NewLeftL+m_NewRightL)*occnum);
            //std::cout<<tempTargetGQN[0].gqn[0]<<std::endl;
        if(i>1) delete m_SuperChain;
        addTwoSites(m_NewLeftL-1,m_NewRightL-1);
        m_SuperChain=new SuperChain(m_NewLeftChain,m_NewRightChain,m_LeftChain,m_RightChain,m_TwoFreeSitesHamiltonian[m_NewLeftL-1]);
        m_SuperChain->TargetGQN=tempTargetGQN;
        //snake::math::printvector(tempTargetGQN);
        m_SuperChain->CalGroundState();
        //if(newleft->base.Dim>m_KeptStatesNum)
        m_SuperChain->renorm(m_KeptStatesNum);
        //else n++;
        m_NewLeftChain->write("./data/L");
        m_NewRightChain->write("./data/R");
        delete m_LeftChain;
        delete m_RightChain;
        m_LeftChain=m_NewLeftChain;
        m_RightChain=m_NewRightChain;
        std::cout<<std::endl;
    }
    delete m_NewLeftChain;
    delete m_NewRightChain;
    delete m_SuperChain;
    //std::cout<<std::endl<<"The maximium exact block size is "<<n<<std::endl<<std::endl;
    std::cout<<"==========Infinite DMRG complete."<<std::endl<<std::endl;
}


void snake::physics::InfiniteDMRG::addTwoSites(int LeftChainLength, int RightChainLength)
{
    m_NewLeftChain=new snake::physics::Chain(*m_LeftChain,snake::physics::allfreesites[LeftChainLength],m_HoppingIntegrals[LeftChainLength-1], m_OnSitePotentials[LeftChainLength], m_TwoSitesInteraction[LeftChainLength-1]);
    std::cout<<"NewLeftL="<<m_NewLeftChain->sitenum<<"\t";
    m_NewRightChain=new snake::physics::Chain(snake::physics::allfreesites[m_ChainLength-RightChainLength-1],
                                       *m_RightChain,m_HoppingIntegrals[m_ChainLength-RightChainLength-1],
                                       m_OnSitePotentials[m_ChainLength-RightChainLength-1],
                                       m_TwoSitesInteraction[m_ChainLength-RightChainLength-1]);
    std::cout<<"NewRightL="<<m_NewRightChain->sitenum<<"\t";
}
