#ifndef GQN_H
#define GQN_H


#include <iostream>
#include <iomanip>
#include <vector>

namespace snake
{

namespace math
{

class GQN
{
 public:
    GQN();
  ~GQN();

          std::vector<int> gqn;


          GQN& operator= (const GQN&)     ;
          GQN& operator+=(const GQN&)     ;
          GQN& operator-=(const GQN& Gvar);
          GQN  operator+ (const GQN&)const;
          GQN  operator- (const GQN&)const;
          bool operator==(const GQN&)const;
          bool operator==(const std::vector<snake::math::GQN>&)const; //Return true if gqn = any gqn in gqnvector. Used for targetting multi gqn number.
          bool operator!=(const GQN&)const;
          bool operator> (const GQN&)const;
          bool operator< (const GQN&)const;

          friend std::ostream& snake::math::operator<<(std::ostream&, const GQN&);
    void none();
    void write(std::ofstream &fout);
    void read(std::ifstream &fin);
    GQN& operator=(int n);
    int num;
};

}
}
#endif
