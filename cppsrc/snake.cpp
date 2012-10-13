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
//ofstream foccnum,;

int main()
{
  double time_start = time(0);

  DMRG chain1;
  chain1.mkdir();

  chain1.iDMRG();

  chain1.fDMRG();
  chain1.iDMRG2tDMRG();
  chain1.tDMRG();

  cout<<endl<<"The multipling times is "<<endl<<multnum<<endl;
  cout << endl << "Total CPU time = " << (time(0) - time_start )  << " seconds" << endl ;
  return EXIT_SUCCESS;
}




