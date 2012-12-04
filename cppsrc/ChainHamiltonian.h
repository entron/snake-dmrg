#ifndef ChainHamiltonian_H
#define ChainHamiltonian_H

/**
Hamiltonian of a block. Sometimes the class can be used as other block operators

@author Cheng Guo
*/


#include "gqnbase.h"
#include "gqnmat.h"
#include "setting.h"

#include<gmd.h>

#include<vector>
#include<fstream>

namespace snake
{

namespace physics
{

class DTMat;
class Site;

class ChainHamiltonian{
public:

/**The block hamiltonian.Make sure that the bases are aligned in the
 *asendent order of the good quantum number.*/
  Rmatrix H;
  Cmatrix HC;

/**If value_type='r', the data are real and use H; if value_type='c', the data *are complex and use HC */
  char value_type;

  ///Base of ChainHamiltonian
  snake::math::GQNBase base;

public:
    ChainHamiltonian();
    ~ChainHamiltonian();

///Block hamiltonian of a block with only one site.Model depended.
  ChainHamiltonian(Site &first,double OnSiteE);

  ///Construct ChainHamiltonian from file "fin"
  ChainHamiltonian(std::ifstream &fin);


  ///Evaluate new ChainHamiltonian
  ChainHamiltonian(ChainHamiltonian *oldham,Site &add,Rmatrix &Hinter,char hand);
  ChainHamiltonian(ChainHamiltonian *oldham,Site &add,Cmatrix &Hinter,char hand);

  ///Evaluate new ChainHamiltonian,useful when calulate local energy
  ChainHamiltonian(ChainHamiltonian *oldham,Site &add,char hand);

  ///Truncate the Block hamiltonian with truncation matrix.
  void renorm(DTMat &mat);

   void write(std::ofstream &fout);
   void read(std::ifstream &fin);

///Convert the data to complex type
  void toComplex();

};

}
}

#endif
