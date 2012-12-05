#ifndef FINITEDMRG_H
#define FINITEDMRG_H

#include "InfiniteDMRG.h"

namespace snake
{
namespace physics
{

class FiniteDMRG: public InfiniteDMRG
{
public:
    FiniteDMRG(InfiniteDMRG &iDMRG);
    void run();
    void CalN();
    snake::physics::SuperChain& generateApdativeTimeDependentDMRGSuperChain();

private:
    void sweep2Right(int StartLeftChainLength, int EndLeftChainLength);
    void sweep2Left(int StartLeftChainLength, int EndLeftChainLength);
    void readSavedLRBlocks(int LeftChainLength);

    int m_sweepTimes;
    int fDMRG_finalNewLeftL;


};
}}
#endif // FINITEDMRG_H
