#ifndef ADAPTIVETIMEDEPENDENTDMRG_H
#define ADAPTIVETIMEDEPENDENTDMRG_H

#include "dmrg.h"
namespace snake
{
namespace physics
{
class AdaptiveTimeDependentDMRG: public DMRG
{
public:
    AdaptiveTimeDependentDMRG(SuperChain &superChain);
    void run();

private:
    void readParameters();

    int m_StepNum;
    //Trotter terms of the real-time evolution operator
    std::vector<LaGenMatComplex> m_LocalUt;
    //The time-dependent Trotter term of the impurity and the first site of the real-time evolution operator
    std::vector<LaGenMatComplex> m_TimeDependentImpurityUt;
};
}}
#endif // ADAPTIVETIMEDEPENDENTDMRG_H
