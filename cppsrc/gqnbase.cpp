#include "gqnbase.h"
#include<algorithm>
#include "public.h"

GQNBase::GQNBase()
{
  Dim=0;
  subnum=0;
}


GQNBase::~GQNBase()
{
}


/*!
\fn GQNBase::kron(GQNBase &b1,GQNBase &b2)
 */

GQNBase kron(const GQNBase& b1,const GQNBase& b2)
{
  GQNBase b;

  b.tempdim.resize(b1.subnum*b2.subnum);
  b.tempsubgqn.resize(b1.subnum*b2.subnum);
  std::vector<GQN> tempsubgqnC(b1.subnum*b2.subnum);

  int n=0;
  for(int i=0;i<b1.subnum;i++)
    for(int j=0;j<b2.subnum;j++)
    {
      b.tempsubgqn[n]=b1.subgqn[i]+b2.subgqn[j];
      tempsubgqnC[n]=b.tempsubgqn[n];
      b.tempdim[n]=b1.dim[i]*b2.dim[j];
      n++;
    }

  std::vector<GQN> subgqnordered(b.tempsubgqn);
  sort(subgqnordered.begin(),subgqnordered.end());
  //for(int i=0;i<subgqnordered.size();i++)
  // std::cout<<subgqnordered[i]<<std::endl;

  ///Caculate map
  b.map.resize(b1.subnum*b2.subnum);
  for(int i=0;i<subgqnordered.size();i++)
    for(int j=0;j<tempsubgqnC.size();j++)
      if(tempsubgqnC[j]==subgqnordered[i])
      {
        b.map[i]=j;
        tempsubgqnC[j].none();///This means pass j in the following step as in ordinary problem good quantum number are not so big
        break;
      }
  // for(int i=0;i<b.map.size();i++) std::cout<<b.map[i]<<" ";  std::cout<<std::endl;

  b.Dim=b1.Dim*b2.Dim;

  ///Caculate subnum
  b.subnum=1;
  GQN lastsubgqn=subgqnordered[0];
  for(int i=1;i<subgqnordered.size();i++)
    if(subgqnordered[i]!=lastsubgqn)
    {
      b.subnum++;
      lastsubgqn=subgqnordered[i];
    }

  b.dim.resize(b.subnum);
  b.subgqn.resize(b.subnum);
  ///Calculate dim and subgqn
  lastsubgqn=subgqnordered[0];
  b.subgqn[0]=lastsubgqn;
  b.dim[0]=b.tempdim[b.map[0]];
  int sub=0;
  for(int i=1;i<subgqnordered.size();i++)
    if(subgqnordered[i]!=lastsubgqn)
    {
      sub++;
      b.subgqn[sub]=subgqnordered[i];
      lastsubgqn=b.subgqn[sub];
      b.dim[sub]=b.tempdim[b.map[i]];
    }
  else
    b.dim[sub]+=b.tempdim[b.map[i]];
  b.genordermap(b1,b2);
  return b;
}


/*!
\fn GQNBase::truncate(GQNBase &dbase,GQNBase &tbase,LaVectorDouble &eigval,double cutedge,int *mark)
 */
void truncate(GQNBase &denmatbase,GQNBase &tmatbase,LaVectorDouble &eigval,double cutedge,int *mark)
{
  int subnum=denmatbase.subnum;
  int m=0,n=0,p=0,q=0;

  //Note that tmatsubnum<=theblock->hamiltonian->subnum
  tmatbase.dim.resize(subnum);
  tmatbase.subgqn.resize(subnum);
  tmatbase.Dim=0;
  for(int i=0;i<subnum;i++)
    tmatbase.dim[i]=0;
  ///Note that some good quantum of the block number may be absent in trunmat.
  for(int i=0;i<subnum;i++)
  {
    tmatbase.subgqn[m]=denmatbase.subgqn[i];
    n=0;
    for(int j=0;j<denmatbase.dim[i];j++)
    {
      if(eigval(p)>cutedge)
      {
        n=1;
        mark[p]=1;///Means the p-th eigenvector are not truncated.
        //trunmat(LaIndex(0,Dim-1),LaIndex(q)).inject(eigvec(LaIndex(0,Dim-1),LaIndex(p)));
        q++;
        tmatbase.dim[m]+=1;
      }
      else mark[p]=0;
      p++;
    }
    tmatbase.Dim+=tmatbase.dim[m];
    if(n==1)m++;
  }
  tmatbase.subnum=m;
}


/*!
\fn GQNBase::write(fstream &fout)
 */
void GQNBase::write(std::ofstream &fout)
{
  fout.write((char*)&Dim,sizeof Dim);
  fout.write((char*)&subnum,sizeof subnum);
  for(int i=0;i<subnum;i++) fout.write((char*)&dim[i],sizeof(int));
  for(int i=0;i<subnum;i++) subgqn[i].write(fout);
//#if USE_INFINITE_DMRG && USE_EVOLVE
  writevector(fout,ordermap);
//#endif
}


/*!
\fn GQNBase::read(std::ifstream &fin)
 */
void GQNBase::read(std::ifstream &fin)
{
  fin.read((char*)&Dim,sizeof Dim);
  fin.read((char*)&subnum,sizeof subnum);
  dim.resize(subnum);
  subgqn.resize(subnum);
  for(int i=0;i<subnum;i++) fin.read((char*)&dim[i],sizeof(int));
  for(int i=0;i<subnum;i++) subgqn[i].read(fin);
//#if USE_INFINITE_DMRG && USE_EVOLVE
  // ordermap.resize(Dim);
  readvector(fin,ordermap);
//#endif
}


/*!
\fn GQNBase::genordermap(GQNBase& b1,GQNBase& b2)
 */
void GQNBase::genordermap(const GQNBase& b1,const GQNBase& b2)
{
  ordermap.resize(Dim);
  int orderrow,start=0;
  for(int k=0;k<subnum;k++)
  {
    orderrow=0;
    for(int l=0;l<b1.subnum;l++)
      for(int m=0;m<b1.dim[l];m++)
        for(int i=0;i<b2.subnum;i++)
          for(int j=0;j<b2.dim[i];j++)
          {
            if(b1.subgqn[l]+b2.subgqn[i]==subgqn[k])
            {
              ordermap[start]=orderrow;
              start++;
            }
            orderrow++;
          }
  }
}

void GQNBase::genvacuumbase()
{
	Dim=1;
	subnum=1;
	dim.resize(1);
	subgqn.resize(1);
	dim[0]=1;
	subgqn[0]=0;

}


/*!
\fn GQNBase::operator=(GQNBase &b)
 */

GQNBase& GQNBase::operator=(const GQNBase& b)
{
  dim=b.dim;
  subgqn=b.subgqn;
  Dim=b.Dim;
  subnum=b.subnum;
  ordermap=b.ordermap;
  tempdim=b.tempdim;
  tempsubgqn=b.tempsubgqn;
  map=b.map;
}

/*!
\fn GQNBase::operator==(GQNBase& base)
 */
int GQNBase::operator==(const GQNBase& base) const
{
  if(dim==base.dim&&Dim==base.Dim&&subgqn==base.subgqn&&subnum==base.subnum)
    return 1;
  else
    return 0;
}


/*!
\fn GQNBase::operator!=(GQNBase& base)
 */
int GQNBase::operator!=(const GQNBase& base) const
{
  if(*this==base)
    return 0;
  else
    return 1;
}

std::ostream & operator<<(std::ostream& os, const GQNBase& base)
{
  std::cout<<"-----GQNBase-----"<<std::endl;
  std::cout<<"Dim="<<base.Dim<<"  subnum="<<base.subnum<<std::endl;
  std::cout<<"subgqn       dim"<<std::endl;
  for(int i=0;i<base.subnum;i++)
  {
    std::cout<<base.subgqn[i]<<"     "<<base.dim[i]<<std::endl;
  }
  return os;
}
