#include "dmrg.h"
#include "dtmat.h"
#include "SuperChain.h"
#include <algorithm>

namespace snake
{
namespace physics
{
    std::vector<snake::physics::Site> allfreesites;
}
}

snake::physics::DMRG::DMRG()
{
    m_KeptStatesNum=50;
    m_TargetGQN.resize(1);
    std::cout<<"==========Parameters:"<<std::endl;
    std::cout<<"KeptStatesNum="<<m_KeptStatesNum<<std::endl;
}


snake::physics::DMRG::~DMRG()
{
		system("rm -rf data");
}



void snake::physics::DMRG::mkdir()
{
	system("rm -rf results");
	system("mkdir results");
	system("rm -rf data");
	system("mkdir data");
}



