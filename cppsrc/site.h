#ifndef SITE_H
#define SITE_H

/**
This class contains the information of the operators of a site.

@author Cheng Guo
*/

#include "gqnmat.h"
#include "setting.h"
#include "gqnbase.h"

#include <gmd.h>
#define LA_COMPLEX_SUPPORT
#include <gmc.h>

#include <vector>
#include <fstream>

namespace snake
{

namespace physics
{


class Block;
class DTMat;

class Site
{
public:
  snake::math::GQNBase base;
  
  ///annilation,creation and number operator.
  std::vector<Rmatrix> a,c,n;
  std::vector<Cmatrix> aC,cC,nC;
  
  ///'r' means real;'c' means complex
  char value_type;

  int num;
  
public:
  Site();
  ~Site();
  ///Read site information from fin
  Site(std::ifstream &fin);
  Site(const Site &s);
  
  ///Evaluate c from a.
  void eval();
  
  ///Renormalize the site operators.
  void renorm(DTMat &mat);
  
  /**Generate the information of sites in block after a new site "add" 
  *is added to the "hand"side of the block.*/
  void addsite(Site &add,char hand);
  
  ///Generate the information of the site after being added to the block.
  void addtoblock(Block &b,char hand);
  
  ///Generate a sigle site.  
 // void genfreesite();
  
  ///Write the site infomation to file  
  void write(std::ofstream &fout);
  void read(std::ifstream &fin);
  //Read site operators from matlab generated files
    void readsite(std::ifstream &basefin, std::ifstream &siteopfin);
	
  Site& operator=(const Site& s);
  
  void toComplex();
    void genfermion();
    void multsignmat();
    void genspin();
  friend std::ostream & operator<<(std::ostream& os, const Site& site);
    void genspinlessfermion();
    void genspinFTDMRG();
};

}
}

#endif

