#ifndef INFINITEDMRG_H
#define INFINITEDMRG_H

#include "dmrg.h"

namespace snake
{
namespace physics
{


class InfiniteDMRG: public DMRG
{
public:
    InfiniteDMRG();
    void run();

protected:
    void addTwoSites(int LeftChainLength, int RightChainLength);
    void readFreesites();
    void readParameters();
    void printParameters();


    Chain *m_LeftChain;
    Chain *m_RightChain;
    Chain *m_NewLeftChain;
    Chain *m_NewRightChain;
    std::vector<double> m_HoppingIntegrals;
    std::vector<double> m_OnSitePotentials;
    std::vector<double> m_TwoSitesInteraction;
    std::vector<LaGenMatDouble> m_TwoFreeSitesHamiltonian;

private:
    void initLeftRightChains();



};


}
}
#endif // INFINITEDMRG_H
