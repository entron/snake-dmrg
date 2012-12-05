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
    void sweep2Right(int StartLeftChainLength, int EndLeftChainLength);
    void sweep2Left(int StartLeftChainLength, int EndLeftChainLength);
    void readSavedLRBlocks(int LeftChainLength);
    void CalN();
    snake::physics::SuperChain& generateApdativeTimeDependentDMRGSuperChain();

    int m_sweepTimes;
private:
    int fDMRG_finalNewLeftL;


};
}}
#endif // FINITEDMRG_H
