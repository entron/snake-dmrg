#ifndef DMRG_H
#define DMRG_H

#include "SuperChain.h"

namespace snake
{

namespace physics
{

class DMRG{
public:
    DMRG();
    ~DMRG();
    void mkdir();
    /// Caculate the average value of onsite operator;
    void CalN();

protected:
    int NewLeftL, NewRightL;
    std::vector<snake::math::GQN> TargetGQN;
    int m_ChainLength;
    int m_KeptStatesNum;
    SuperChain *m_SuperChain;
};

}
}
#endif
