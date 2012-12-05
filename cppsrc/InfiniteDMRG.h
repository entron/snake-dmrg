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
    void readFreesites();
    void readParameters();
    void run();
    void addTwoSites(int LeftChainLength, int RightChainLength);


protected:
    void ReadSavedLRBlocks(int LeftChainLength);
    void AddTwoSites(int LeftChainLength, int RightChainLength);
    Chain *left,*right,*newleft,*newright;
    std::vector<double> HoppingIntegrals;
    std::vector<double> OnSitePotentials;
    std::vector<double> TwoSitesInteraction;
    std::vector<LaGenMatDouble> H;

};


}
}
#endif // INFINITEDMRG_H
