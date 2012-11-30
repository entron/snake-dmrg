#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

using namespace std;
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


long int multnum=0; //Using for debug
//Site freesite;
double totaltrunerror=0;
vector<Site> allfreesites;
//std::ofstream foccnum,;

int main()
{
  double time_start = time(0);

  DMRG chain1;
  chain1.mkdir();

  chain1.iDMRG();

  chain1.fDMRG();
  chain1.iDMRG2tDMRG();
  chain1.tDMRG();

  std::cout<<std::endl<<"The multipling times is "<<std::endl<<multnum<<std::endl;
  std::cout << std::endl << "Total CPU time = " << (time(0) - time_start )  << " seconds" << std::endl ;
  return EXIT_SUCCESS;
}




