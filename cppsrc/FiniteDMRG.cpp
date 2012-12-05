#include "FiniteDMRG.h"

snake::physics::FiniteDMRG::FiniteDMRG(
        snake::physics::InfiniteDMRG& iDMRG): InfiniteDMRG(iDMRG)
{
    m_sweepTimes = 1;
    std::cout<<"fDMRGSweepTimes="<<m_sweepTimes<<std::endl;
}

void snake::physics::FiniteDMRG::run()
{
    std::cout<<"==========Start Finite DMRG."<<std::endl;
    fDMRG_finalNewLeftL=m_NewLeftL;
    for(int num=0;num<m_sweepTimes;num++)
    {
        //left increase at the cost of right
        std::cout<<"*****No."<<num<<" sweep"<<std::endl;
        sweep2Left(m_NewLeftL-1,1);
        sweep2Right(1,m_ChainLength-3);
        sweep2Left(m_ChainLength-3,fDMRG_finalNewLeftL-1);
    }
    std::cout<<"==========Finite DMRG complete."<<std::endl<<std::endl;
}

void snake::physics::FiniteDMRG::sweep2Right(int StartLeftChainLength, int EndLeftChainLength)
{
    std::cout<<"-----Start sweeping to right."<<std::endl;
    for(int i=StartLeftChainLength;i<=EndLeftChainLength;i++)
    {
        readSavedLRBlocks(i);
        addTwoSites(i,m_ChainLength-i-2);
        m_SuperChain=new SuperChain(m_NewLeftChain,m_NewRightChain,m_LeftChain,m_RightChain,m_TwoFreeSitesHamiltonian[i]);
        m_SuperChain->TargetGQN=m_TargetGQN;
        m_SuperChain->CalGroundState();
        //	if(newleft->base.Dim>m_KeptStatesNum)
        m_SuperChain->renormleft(m_KeptStatesNum);
        m_NewLeftChain->write("./data/L");
        delete m_LeftChain;
        delete m_RightChain;
        delete m_NewLeftChain;
        delete m_NewRightChain;
        delete m_SuperChain;
        std::cout<<std::endl;
    }
}


void snake::physics::FiniteDMRG::sweep2Left(int StartLeftChainLength, int EndLeftChainLength)
{
    std::cout<<"-----Start sweeping to left."<<std::endl;
    for(int i=StartLeftChainLength;i>=EndLeftChainLength;i--)
    {
        readSavedLRBlocks(i);
        addTwoSites(i,m_ChainLength-i-2);
        m_SuperChain=new SuperChain(m_NewLeftChain,m_NewRightChain,m_LeftChain,m_RightChain,m_TwoFreeSitesHamiltonian[i]);
        m_SuperChain->TargetGQN=m_TargetGQN;
        m_SuperChain->CalGroundState();
        //std::cout<<m_SuperChain->wf<<std::endl;break;
        //	if(newright->base.Dim>m_KeptStatesNum)
        m_SuperChain->renormright(m_KeptStatesNum);

        m_NewRightChain->write("./data/R");
        delete m_LeftChain;
        delete m_RightChain;
        delete m_NewRightChain;
        delete m_NewLeftChain;
        delete m_SuperChain;
        std::cout<<std::endl;
    }
}

snake::physics::SuperChain& snake::physics::FiniteDMRG::generateApdativeTimeDependentDMRGSuperChain()
{
    std::cout<<"==========Passing fDMRG m_SuperChain to the t-DMRG."<<std::endl;
    ///Pass the DMRG m_SuperChain to the tDMRG
    readSavedLRBlocks(fDMRG_finalNewLeftL-1);
    addTwoSites(fDMRG_finalNewLeftL-1,fDMRG_finalNewLeftL-2);
    m_SuperChain=new SuperChain(m_NewLeftChain,m_NewRightChain,m_LeftChain,m_RightChain,m_TwoFreeSitesHamiltonian[fDMRG_finalNewLeftL-1]);
    //TargetGQN.resize(2);
    m_SuperChain->TargetGQN=m_TargetGQN;
    m_SuperChain->CalGroundState();
    if(m_NewLeftChain->base.Dim>m_KeptStatesNum)
        m_SuperChain->renorm(m_KeptStatesNum);
    m_SuperChain->KeptStatesNum=m_KeptStatesNum;
    //std::cout<<freesite<<std::endl;
    m_SuperChain->loaddtmats();
    m_SuperChain->sweep2left(fDMRG_finalNewLeftL);
    //////////////////////////////
    m_SuperChain->TargetGQN2.resize(1);
    m_SuperChain->TargetGQN2[0].gqn[0]=m_SuperChain->TargetGQN[0].gqn[0];


    ///Generate new wave functions
    Rmatrix wfmat2=m_SuperChain->wfmat;
    int row=wfmat2.rowbase.subnum;
    int col=wfmat2.colbase.subnum;
    for(int i=0;i<row;i++)
        if(wfmat2.pmat(i,0)!=-1)
        {
            wfmat2.pmat(i,1)=wfmat2.pmat(i,0);
            wfmat2.pmat(i,0)=-1;
        }
    m_SuperChain->extractwf(m_SuperChain->wfmat,m_SuperChain->wf,m_SuperChain->TargetGQN);
    m_SuperChain->toComplex();
    R2C(wfmat2,m_SuperChain->wfmatC2);  m_SuperChain->extractwf(m_SuperChain->wfmatC2,m_SuperChain->wfC2,m_SuperChain->TargetGQN2);
    //std::cout<<wfmat1<<std::endl;std::cout<<wfmat2<<std::endl;
    for(int i=0;i<m_ChainLength;i++)
    {
        snake::physics::allfreesites[i].toComplex();
        snake::physics::allfreesites[i].eval();
    }
    delete m_LeftChain;
    delete m_RightChain;
    delete m_NewLeftChain;
    delete m_NewRightChain;
    std::cout<<std::endl;
    system("rm -rf data");
    std::cout<<"==========Passing fDMRG m_SuperChain to the t-DMRG complete."<<std::endl<<std::endl;
    return *m_SuperChain;
}


void snake::physics::FiniteDMRG::readSavedLRBlocks(int LeftChainLength)
{
    std::string fname;
    std::stringstream stl,str;
    fname="./data/L";
    stl<<LeftChainLength;
    fname=fname+stl.str();
    // std::cout<<fname<<std::endl;
    m_LeftChain=new snake::physics::Chain(fname);
    fname="./data/R";
    str<<(m_ChainLength-LeftChainLength-2);
    fname=fname+str.str();
    //std::cout<<fname<<std::endl;
    m_RightChain=new snake::physics::Chain(fname);
}


void snake::physics::FiniteDMRG::CalN()
{
    std::cout<<"Calculate N of Everysites when Chain sweeped to the middle."<<std::endl;
    readSavedLRBlocks(fDMRG_finalNewLeftL-1);
    addTwoSites(fDMRG_finalNewLeftL-1,fDMRG_finalNewLeftL-2);
    m_SuperChain=new SuperChain(m_NewLeftChain,m_NewRightChain,m_LeftChain,m_RightChain,m_TwoFreeSitesHamiltonian[fDMRG_finalNewLeftL-1]);
    m_SuperChain->TargetGQN=m_TargetGQN;
    m_SuperChain->CalGroundState();
    if(m_NewLeftChain->base.Dim>m_KeptStatesNum)
        m_SuperChain->renorm(m_KeptStatesNum);

    //m_SuperChain->calN("./data/n.dat");
    std::ofstream fout("./data/n.dat");
    m_NewLeftChain->calN(fout,'l');
    m_NewRightChain->calN(fout,'r');
    fout.close();

    system("xmgrace ./data/n.dat &");
    //m_SuperChain->prepare(); ///Store base information of the blocks
    delete m_LeftChain;
    delete m_RightChain;
    delete m_NewLeftChain;
    delete m_NewRightChain;
    delete m_SuperChain;
}
