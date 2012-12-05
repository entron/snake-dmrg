#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

#include<gmd.h>
#include<lavd.h>
#include"site.h"
#include"Chain.h"
#include"SuperChain.h"
#include "ChainHamiltonian.h"
#include "public.h"
#include "setting.h"
#include <blaspp.h>
#include "dtmat.h"
#include "gqnmat.h"
#include <string>
#include <sstream>
#include "InfiniteDMRG.h"
#include "FiniteDMRG.h"
#include "AdaptiveTimeDependentDMRG.h"

namespace snake
{
long int multnum=0; //Using for debug
//std::vector<snake::physics::Site> allfreesites;
}

//Site freesite;

//std::ofstream foccnum,;


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

  std::cout<<std::endl<<"The multipling times is "<<std::endl<<snake::multnum<<std::endl;
  std::cout << std::endl << "Total CPU time = " << (time(0) - time_start )  << " seconds" << std::endl ;
  return EXIT_SUCCESS;
}




