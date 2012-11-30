#ifndef BLOCK_H
#define BLOCK_H

/**
This class contains information related to a block.

@author Cheng Guo
*/

#include "gqnbase.h"
#include "setting.h"
#include "gqnmat.h"

#include<gmd.h>

#include<vector>
#include<fstream>
#include<string>

class Site;
class BlocHam;
class DTMat;

class Block{
public:
  int sitenum;

  ///Information of each site in the block;
  std::vector<Site> site;

  ///For the convenice of SupHam
  Site *siteadded;

  BlocHam *hamiltonian;

  DTMat *dtmat;

  GQNBase base;

  ///The interaction term of the right(left) most site and the new added site.
  Rmatrix Hinter;
  Cmatrix HinterC;

  char value_type;

public:

  Block();
  ~Block();

  ///Generate block with only one site
  Block(Site &first,double OnSiteE);

  ///New block with the site added to the right of the old block.
  Block(Block &old,Site &add,double HoppingT, double OnSiteE, double TwoSitesV);

  ///New block with the site added to the left of the old block.
  //Block(Site &add,Block &old);
	Block(Site &addsite,Block &old,double HoppingT, double OnSiteE, double TwoSitesV);

   /**The block hamiltonian is only the local energy of localsite
  *and localsite+1.This function is used in time evolution*/
  //Block(Site& addsite,Block& old,int localsite);

  ///Construct Block from file
  Block(std::string &filename);

  /**This function calculate the interaction hamiiltonian of two site,
   *and is model depended.*/
  //void calHinter(Site &a,Site &b,char addsiteposition,char include_onsite='y');
	void calHinter(Site &siteA,Site &siteB,char addsiteposition,double HoppingT, double OnSiteE, double TwoSitesV);

  void renorm();

  ///Write the block information to a file,whose name starts with prefix
  void write(char* prefix);
  void read(std::string &filename);

  ///Calculat the average value of operator n at each site.
  void calN(std::ofstream &fout,char hand);

  void toComplex();
    void calHinter_Heisenberg(Site &siteA,Site &siteB,char addsiteposition);

private:

  void initialadd(Block &old,Site &add);
  ///Renorm all sites
  void renormsites();
  ///Renorm only the side site
  void renormsidesite();
  void calHinter_Hubbard(Site &siteA,Site &siteB,char addsiteposition,char include_onsite='y');

};


#endif
