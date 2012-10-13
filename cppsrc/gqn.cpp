#include "setting.h"
#include "gqn.h"

GQN::GQN()
{
  num=NUMBER_OF_GOOD_QUANTUM_NUMBER;
  gqn.resize(num);
}

GQN::~GQN()
{

}


GQN& GQN::operator= (const GQN& Gvar)
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       gqn[i]=Gvar.gqn[i];
    }
    return *this;
}
bool GQN::operator==(const GQN& Gvar)const
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       if(gqn[i]!=Gvar.gqn[i]) return false;
    }
    return true;
}
bool GQN::operator==(const vector<GQN>& Gvecvar)const
{
	int vlen=Gvecvar.size();
	for(int n=0;n<vlen;n++)
	{
		if(gqn==Gvecvar[n].gqn) return true;
	}
    return false;
}


bool GQN::operator!=(const GQN& Gvar)const
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       if(gqn[i]!=Gvar.gqn[i]) return true;
    }
    return false;
}
bool GQN::operator> (const GQN& Gvar)const
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       if(gqn[i]>Gvar.gqn[i])return true ;
       if(gqn[i]<Gvar.gqn[i])return false;
    }
    return false;
}
bool GQN::operator< (const GQN& Gvar)const
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       if(gqn[i]<Gvar.gqn[i])return true;
       if(gqn[i]>Gvar.gqn[i])return false;
    }
    return false;
}
GQN& GQN::operator+=(const GQN& Gvar)
{
    for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
    {
       gqn[i]+=Gvar.gqn[i];
    }
    return *this;
}

GQN& GQN::operator-=(const GQN& Gvar)
{
  for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER;i++)
  {
    gqn[i]-=Gvar.gqn[i];
  }
  return *this;
}

GQN  GQN::operator+ (const GQN& Gvar)const
{
    GQN Gsum=*this;
    Gsum+=Gvar;
    return Gsum;
}


GQN  GQN::operator- (const GQN& Gvar)const
{
  GQN Gminus=*this;
  Gminus-=Gvar;
  return Gminus;
}

ostream& operator<<(ostream& output,const GQN& GQvar)
{
    output<<"(";
    if(NUMBER_OF_GOOD_QUANTUM_NUMBER==1)
    {
       output<<setw(3)<<GQvar.gqn[NUMBER_OF_GOOD_QUANTUM_NUMBER-1]<<")   ";
    }
    else
    {
       for(int i=0;i<NUMBER_OF_GOOD_QUANTUM_NUMBER-1;i++)
       {
          output<<setw(3)<<GQvar.gqn[i]<<",";
       }
       output<<setw(3)<<GQvar.gqn[NUMBER_OF_GOOD_QUANTUM_NUMBER-1]<<"  )";
    }
    return output;
}

/*!
    \fn GQN::none()
 */
void GQN::none()
{
  for(int i=0;i<num;i++)
    gqn[i]=99999;
}


/*!
    \fn GQN::write(ofstream &fout)
 */
void GQN::write(ofstream &fout)
{
  fout.write((char*)&num,sizeof num);
  for(int i=0;i<num;i++)
    fout.write((char*)&gqn[i],sizeof gqn[i]);
}


/*!
    \fn GQN::read(ifstream &fin)
 */
void GQN::read(ifstream &fin)
{
  fin.read((char*)&num,sizeof num);
  for(int i=0;i<num;i++)
    fin.read((char*)&gqn[i],sizeof gqn[i]);
}


/*!
    \fn GQN::operator=(int n)
///For the compatibility with earlier codes
 */
GQN& GQN::operator=(int n)
{
  if(num==1)
    gqn[0]=n;
  else
    cout<<"There are more than one good quantum number!"<<endl;
}
