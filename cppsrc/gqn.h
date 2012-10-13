#ifndef GQN_H
#define GQN_H


using namespace std;

#include <iostream>
#include <iomanip>
#include <vector>

class GQN
{
 public:
    GQN();
  ~GQN();

          vector<int> gqn;


          GQN& operator= (const GQN&)     ;
          GQN& operator+=(const GQN&)     ;
          GQN& operator-=(const GQN& Gvar);
          GQN  operator+ (const GQN&)const;
          GQN  operator- (const GQN&)const;
          bool operator==(const GQN&)const;
          bool operator==(const vector<GQN>&)const; //Return true if gqn = any gqn in gqnvector. Used for targetting multi gqn number.
          bool operator!=(const GQN&)const;
          bool operator> (const GQN&)const;
          bool operator< (const GQN&)const;

          friend ostream& operator<<(ostream&, const GQN&);
    void none();
    void write(ofstream &fout);
    void read(ifstream &fin);
    GQN& operator=(int n);
    int num;
};
#endif
