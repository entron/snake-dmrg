#ifndef GQNBASE_H
#define GQNBASE_H

#include <vector>
#include <fstream>
#include <lavd.h>
#include "gqn.h"

/**
This class stores and process the information about good quantum number bases.

@author Guo Cheng
*/

namespace snake
{

namespace math
{

class GQNBase{
public:
 ///The total dim of the base
  int Dim;

  ///Subspace number
  int subnum;

  ///Dimension of each subspace
  std::vector<int> dim;

  ///Good quantum number corresponding to each subspace
  std::vector<snake::math::GQN> subgqn;

  ///Useful when the object of this class is generate by a kronecker product.
  std::vector<int> map;

    /**ordermap is used to rearrange the order of newham's base vectors so
  *that they are ascendent according to newham->TargetGQN.
  */
  std::vector<int> ordermap;

  std::vector<int> tempdim;
  std::vector<snake::math::GQN> tempsubgqn;

public:
  GQNBase();
  ~GQNBase();

  void write(std::ofstream &fout);
  void read(std::ifstream &fin);

  void genordermap(const GQNBase& b1,const GQNBase& b2);
  void genvacuumbase();

  GQNBase& operator=(const GQNBase& b);
  int operator==(const GQNBase& base) const;
  int operator!=(const GQNBase& base) const;

  friend GQNBase kron(const GQNBase &b1,const GQNBase &b2);
  friend void truncate(GQNBase &dbase,GQNBase &tbase,LaVectorDouble &eigval,double cutedge,int *mark);
  friend std::ostream & snake::math::operator<<(std::ostream& os, const GQNBase& base);

};

   GQNBase kron(const GQNBase &b1,const GQNBase &b2);
   void truncate(GQNBase &dbase,GQNBase &tbase,LaVectorDouble &eigval,double cutedge,int *mark);

}
}

#endif
