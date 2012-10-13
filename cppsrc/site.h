#ifndef SITE_H
#define SITE_H

/**
This class contains the information of the operators of a site.

@author Cheng Guo
*/

using namespace std;

#include "gqnmat.h"
#include "setting.h"
#include "gqnbase.h"

#include <gmd.h>
#define LA_COMPLEX_SUPPORT
#include <gmc.h>

#include <vector>
#include <fstream>

class Block;
class DTMat;

class Site
{
public:
  GQNBase base;
  
  ///annilation,creation and number operator.
  vector<Rmatrix> a,c,n; 
  vector<Cmatrix> aC,cC,nC; 
  
  ///'r' means real;'c' means complex
  char value_type;

  int num;
  
public:
  Site();
  ~Site();
  ///Read site information from fin
  Site(ifstream &fin);
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
  void write(ofstream &fout);
  void read(ifstream &fin);
  //Read site operators from matlab generated files
	void readsite(ifstream &basefin, ifstream &siteopfin);
	
  Site& operator=(const Site& s);
  
  void toComplex();
    void genfermion();
    void multsignmat();
    void genspin();
  friend ostream & operator<<(ostream& os, const Site& site);
    void genspinlessfermion();
    void genspinFTDMRG();
};

#endif

