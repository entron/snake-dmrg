#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

#include<gmd.h>
#include<lavd.h>
#include"site.h"
#include"block.h"
#include"supblock.h"
#include "blocham.h"
#include "public.h"
#include "setting.h"
#include <blaspp.h>
#include "dtmat.h"
#include "gqnmat.h"
#include <string>
#include <sstream>
#include "dmrg.h"


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

  snake::physics::DMRG chain1;
  chain1.mkdir();

  chain1.iDMRG();

  chain1.fDMRG();
  chain1.iDMRG2tDMRG();
  chain1.tDMRG();

  std::cout<<std::endl<<"The multipling times is "<<std::endl<<snake::multnum<<std::endl;
  std::cout << std::endl << "Total CPU time = " << (time(0) - time_start )  << " seconds" << std::endl ;
  return EXIT_SUCCESS;
}




