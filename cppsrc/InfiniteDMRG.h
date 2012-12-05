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

    Chain *left,*right,*newleft,*newright;
    std::vector<double> HoppingIntegrals;
    std::vector<double> OnSitePotentials;
    std::vector<double> TwoSitesInteraction;
    std::vector<LaGenMatDouble> H;

};


}
}
#endif // INFINITEDMRG_H
