#ifndef GQNBASE_H
#define GQNBASE_H

#include <vector>
#include <fstream>
using namespace std;
#include <lavd.h>
#include "gqn.h"

/**
This class stores and process the information about good quantum number bases.

@author Guo Cheng
*/

class GQNBase{
public:
 ///The total dim of the base
  int Dim;

  ///Subspace number
  int subnum;

  ///Dimension of each subspace
  vector<int> dim;

  ///Good quantum number corresponding to each subspace
  vector<GQN> subgqn;

  ///Useful when the object of this class is generate by a kronecker product.
  vector<int> map;

    /**ordermap is used to rearrange the order of newham's base vectors so
  *that they are ascendent according to newham->TargetGQN.
  */
  vector<int> ordermap;

  vector<int> tempdim;
  vector<GQN> tempsubgqn;

public:
  GQNBase();
  ~GQNBase();

  void write(ofstream &fout);
  void read(ifstream &fin);

  void genordermap(const GQNBase& b1,const GQNBase& b2);
  void genvacuumbase();

  GQNBase& operator=(const GQNBase& b);
  int operator==(const GQNBase& base) const;
  int operator!=(const GQNBase& base) const;

  friend GQNBase kron(const GQNBase &b1,const GQNBase &b2);
  friend void truncate(GQNBase &dbase,GQNBase &tbase,LaVectorDouble &eigval,double cutedge,int *mark);
  friend ostream & operator<<(ostream& os, const GQNBase& base);

};

#endif
