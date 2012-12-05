#include "AdaptiveTimeDependentDMRG.h"

snake::physics::AdaptiveTimeDependentDMRG::AdaptiveTimeDependentDMRG(snake::physics::SuperChain& superChain)
{
    m_ChainLength = superChain.sitenum;
    m_SuperChain = &superChain;
}

void snake::physics::AdaptiveTimeDependentDMRG::readParameters()
{
    std::cout<<"==========Reading model parameters for t-DMRG from Matlab generated files."<<std::endl;
    std::string fname;
    fname="./model/rt_T0.dat";
    std::ifstream fin(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    fin.read((char*)&m_StepNum,sizeof(int));
    snake::math::ReadOperators(fin,rt_OP,m_ChainLength-1);
    fin.close();
    fname="./model/rt_H1_T0.dat";
    fin.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    snake::math::ReadOperators(fin,rt_td_impurity_OP,m_StepNum);
    fin.close();

    std::cout<<"tDMRGStepNum="<<m_StepNum<<std::endl;
}


void snake::physics::AdaptiveTimeDependentDMRG::run()
{
    readParameters();
    std::cout<<"==========Start t-DMRG."<<std::endl;
    m_SuperChain->creatoutputfiles();
    m_SuperChain->evolve(rt_td_impurity_OP,rt_OP,m_StepNum);
    m_SuperChain->closeoutputfiles();
    delete m_SuperChain;
    std::cout<<"==========t-DMRG complete."<<std::endl<<std::endl;
}
