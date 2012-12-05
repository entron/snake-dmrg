#include <iostream>

#include "InfiniteDMRG.h"
#include "FiniteDMRG.h"
#include "AdaptiveTimeDependentDMRG.h"


int main()
{
  double time_start = time(0);

  snake::physics::InfiniteDMRG iDMRG;
  iDMRG.mkdir();
  iDMRG.run();

  snake::physics::FiniteDMRG fDMRG(iDMRG);
  fDMRG.run();

  snake::physics::AdaptiveTimeDependentDMRG tDMRG(fDMRG.generateApdativeTimeDependentDMRGSuperChain());
  tDMRG.run();

  std::cout << std::endl << "Total CPU time = " << (time(0) - time_start )  << " seconds" << std::endl ;
  return EXIT_SUCCESS;
}




