#include "setting.h"
#include "gqn.h"

snake::math::GQN::GQN()
{
  num=NUMBER_OF_GOOD_QUANTUM_NUMBER;
  gqn.resize(num);
}

snake::math::GQN::~GQN()
{

}


snake::math::GQN& snake::math::GQN::operator= (const snake::math::GQN& Gvar)
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       gqn[i]=Gvar.gqn[i];
    }
    return *this;
}
bool snake::math::GQN::operator==(const snake::math::GQN& Gvar)const
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       if(gqn[i]!=Gvar.gqn[i]) return false;
    }
    return true;
}
bool snake::math::GQN::operator==(const std::vector<snake::math::GQN>& Gvecvar)const
{
	int vlen=Gvecvar.size();
	for(int n=0;n<vlen;n++)
	{
		if(gqn==Gvecvar[n].gqn) return true;
	}
    return false;
}


bool snake::math::GQN::operator!=(const snake::math::GQN& Gvar)const
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       if(gqn[i]!=Gvar.gqn[i]) return true;
    }
    return false;
}
bool snake::math::GQN::operator> (const snake::math::GQN& Gvar)const
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       if(gqn[i]>Gvar.gqn[i])return true ;
       if(gqn[i]<Gvar.gqn[i])return false;
    }
    return false;
}
bool snake::math::GQN::operator< (const snake::math::GQN& Gvar)const
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       if(gqn[i]<Gvar.gqn[i])return true;
       if(gqn[i]>Gvar.gqn[i])return false;
    }
    return false;
}
snake::math::GQN& snake::math::GQN::operator+=(const snake::math::GQN& Gvar)
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       gqn[i]+=Gvar.gqn[i];
    }
    return *this;
}

snake::math::GQN& snake::math::GQN::operator-=(const snake::math::GQN& Gvar)
{
  for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
  {
    gqn[i]-=Gvar.gqn[i];
  }
  return *this;
}

snake::math::GQN  snake::math::GQN::operator+ (const snake::math::GQN& Gvar)const
{
    snake::math::GQN Gsum=*this;
    Gsum+=Gvar;
    return Gsum;
}


snake::math::GQN  snake::math::GQN::operator- (const snake::math::GQN& Gvar)const
{
  snake::math::GQN Gminus=*this;
  Gminus-=Gvar;
  return Gminus;
}

std::ostream& snake::math::operator<<(std::ostream& output,const snake::math::GQN& GQvar)
{
    output<<"(";
    if(NUMBER_OF_GOOD_QUANTUM_NUMBER==1)
    {
       output<<std::setw(3)<<GQvar.gqn[NUMBER_OF_GOOD_QUANTUM_NUMBER-1]<<")   ";
    }
    else
    {
       for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER-1;i++)
       {
          output<<std::setw(3)<<GQvar.gqn[i]<<",";
       }
       output<<std::setw(3)<<GQvar.gqn[NUMBER_OF_GOOD_QUANTUM_NUMBER-1]<<"  )";
    }
    return output;
}


void snake::math::GQN::none()
{
  for(int i=0;i<num;i++)
    gqn[i]=99999;
}



void snake::math::GQN::write(std::ofstream &fout)
{
  fout.write((char*)&num,sizeof num);
  for(int i=0;i<num;i++)
    fout.write((char*)&gqn[i],sizeof gqn[i]);
}


void snake::math::GQN::read(std::ifstream &fin)
{
  fin.read((char*)&num,sizeof num);
  for(int i=0;i<num;i++)
    fin.read((char*)&gqn[i],sizeof gqn[i]);
}



snake::math::GQN& snake::math::GQN::operator=(int n)
{
  if(num==1)
    gqn[0]=n;
  else
    std::cout<<"There are more than one good quantum number!"<<std::endl;
}
