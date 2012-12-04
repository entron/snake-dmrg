#ifndef BLOCHAM_H
#define BLOCHAM_H

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

class BlocHam{
public:

/**The block hamiltonian.Make sure that the bases are aligned in the
 *asendent order of the good quantum number.*/
  Rmatrix H;
  Cmatrix HC;

/**If value_type='r', the data are real and use H; if value_type='c', the data *are complex and use HC */
  char value_type;

  ///Base of BlocHam
  snake::math::GQNBase base;

public:
    BlocHam();
    ~BlocHam();

///Block hamiltonian of a block with only one site.Model depended.
  BlocHam(Site &first,double OnSiteE);

  ///Construct BlocHam from file "fin"
  BlocHam(std::ifstream &fin);


  ///Evaluate new BlocHam
  BlocHam(BlocHam *oldham,Site &add,Rmatrix &Hinter,char hand);
  BlocHam(BlocHam *oldham,Site &add,Cmatrix &Hinter,char hand);

  ///Evaluate new BlocHam,useful when calulate local energy
  BlocHam(BlocHam *oldham,Site &add,char hand);

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
