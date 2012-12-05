#ifndef DMRG_H
#define DMRG_H

#include "Chain.h"
#include "SuperChain.h"
#include "setting.h"
#include "gqn.h"
#include "site.h"


/**
The frontend of SNAKE.

	@author Cheng.Guo <Cheng.Guo@physik.lmu.de>

*/
namespace snake
{

namespace physics
{

class DMRG{
public:
    DMRG();

    ~DMRG();
    void mkdir();
    ///Caculate the average value of onsite operator;
    void CalN();

    int m_ChainLength;
    int m_KeptStatesNum;
    SuperChain *m_SuperChain;


protected:
	int NewLeftL, NewRightL;
  ///The number of site which are exact.
  //int n;
  std::vector<snake::math::GQN> TargetGQN;

};

}
}
#endif
