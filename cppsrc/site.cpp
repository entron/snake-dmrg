#include "site.h"
#include "block.h"
#include "public.h"
#include "dtmat.h"
#include "setting.h"

#include <iostream>

#include <blas3pp.h>
#include <blaspp.h>


Site::Site()
{
  value_type='r';
  num=NUMBER_OF_KINDS_OF_PARTICLES;

}

Site::~Site()
{
}

/*!
\fn Site::Site(const Site &s)
 */
Site::Site(const Site &s)
{
    /// @todo implement me
}

void Site::eval()
{
  if(value_type=='r')
  {
    Rmatrix identity;
    c.resize(num);
    for(int i=0;i<num;i++)
    {
      c[i].resize(base,base);
      identity.geneye(a[i]);
      Mat_Trans_Mat_Mult(a[i],identity,c[i]);
    }
  }
  else
  {
    Cmatrix identity;
    cC.resize(num);
    for(int i=0;i<num;i++)
    {
      cC[i].resize(base,base);
      identity.geneye(aC[i]);
      Mat_Trans_Mat_Mult(aC[i],identity,cC[i]);
    }
}
}

/*!
    \fn Site::renorm(LaGenMat &trunmat)
 */
void Site::renorm(DTMat &mat)
{
  int KeptStatesNum;
  KeptStatesNum=mat.tmatbase.Dim;
  if(value_type=='r')
  {
  Rmatrix midx(mat.tmatbase,base),midy(mat.tmatbase,mat.tmatbase);
    for(int i=0;i<num;i++)
    {
      Mat_Trans_Mat_Mult(mat.trunmat,a[i],midx);
      Mat_Mat_Mult(midx,mat.trunmat,midy);
      a[i]=midy;
      Mat_Trans_Mat_Mult(mat.trunmat,n[i],midx);
      Mat_Mat_Mult(midx,mat.trunmat,midy);
      n[i]=midy;
    }
  }
  else
  {
    Cmatrix midx(mat.tmatbase,base),midy(mat.tmatbase,mat.tmatbase);
    //std::cout<<mat.trunmattrans2<<std::endl;
    //std::cout<<aC<<std::endl;
    for(int i=0;i<num;i++)
    {
      Mat_Trans_Mat_Mult(mat.trunmatC,aC[i],midx);
      Mat_Mat_Mult(midx,mat.trunmatC,midy);
      aC[i]=midy;
      Mat_Trans_Mat_Mult(mat.trunmatC,nC[i],midx);
      Mat_Mat_Mult(midx,mat.trunmatC,midy);
      nC[i]=midy;
    }
  }
  base=mat.tmatbase;
}


/*!
    \fn Site::addsite(Site & add,char hand)
 */
void Site::addsite(Site & add,char hand)
{
  if(value_type=='r')
  {
    Rmatrix identity,temp;
    identity.geneye(add.base);
  //add.geniden();
    if(hand=='r')
    {
      for(int i=0;i<num;i++)
      {
        kron(a[i],identity,temp);
        a[i]=temp;
        kron(n[i],identity,temp);
        n[i]=temp;
      }
      base=kron(base,add.base);
    }
    else
    {
      for(int i=0;i<num;i++)
      {
        kron(identity,a[i],temp);
        a[i]=temp;
        kron(identity,n[i],temp);
        n[i]=temp;
      }
      base=kron(add.base,base);
    }

  }
  else
  {
    Cmatrix identity,temp;
    identity.geneye(add.base);
  //add.geniden();
    if(hand=='r')
    {
      for(int i=0;i<num;i++)
      {
        kron(aC[i],identity,temp);
        aC[i]=temp;
        kron(nC[i],identity,temp);
        nC[i]=temp;
      }
      base=kron(base,add.base);
    }
    else
    {
      for(int i=0;i<num;i++)
      {
        kron(identity,aC[i],temp);
        aC[i]=temp;
        kron(identity,nC[i],temp);
        nC[i]=temp;
      }
      base=kron(add.base,base);
    }

  }
}



/*!
    \fn Site::addtoblock(Block &b,char hand)
 */
void Site::addtoblock(Block &b,char hand)
{
  if(value_type=='r')
  {
    Rmatrix identity,temp;
    identity.geneye(b.base);

    if(hand=='r')
    {
#if FERMIONSIGN
      for(int i=0;i<num;i++)
      {
        Rmatrix signmat,tempmat;
        signmat.gensignmat(base,i);
        tempmat=a[i]*signmat;
        //std::cout<<identity<<std::endl;
        //std::cout<<tempmat<<std::endl;
        kron(identity,tempmat,a[i]);
        kron(identity,n[i],temp);
        n[i]=temp;
      }
#else
      for(int i=0;i<num;i++)
      {
        kron(identity,a[i],temp);
        a[i]=temp;
        kron(identity,n[i],temp);
        n[i]=temp;
      }
#endif
      base=kron(b.base,base);
    }
    else
    {
      for(int i=0;i<num;i++)
      {
        kron(a[i],identity,temp);
        a[i]=temp;
        kron(n[i],identity,temp);
        n[i]=temp;
      }
      base=kron(base,b.base);
    }
  }
  else
  {
    Cmatrix identity,temp;
    identity.geneye(b.base);

    if(hand=='r')
    {
      for(int i=0;i<num;i++)
      {
        kron(identity,aC[i],temp);
        aC[i]=temp;
        kron(identity,nC[i],temp);
        nC[i]=temp;
      }
      base=kron(b.base,base);
    }
    else
    {
      for(int i=0;i<num;i++)
      {
        kron(aC[i],identity,temp);
        aC[i]=temp;
        kron(nC[i],identity,temp);
        nC[i]=temp;
      }
      base=kron(base,b.base);
    }
  }

}


/*
void Site::genfreesite()
{

	std::string fname="./model/site_operators.dat";
	std::ifstream opfin(fname.c_str(),std::ios_base::in|std::ios_base::binary);
	std::string fname2="./model/site_base.dat";
	std::ifstream basefin(fname2.c_str(),std::ios_base::in|std::ios_base::binary);

	readsite(basefin, opfin);
	opfin.close();
	basefin.close();
   // genspinFTDMRG();
  //genspin();
  //genspinlessfermion();
  //genfermion();
	//std::cout<<"==========Free site information:"<<std::endl;
	//std::cout<<*this;
}
*/

void Site::readsite(std::ifstream &basefin, std::ifstream &siteopfin)
{
  value_type='r';
  num=NUMBER_OF_KINDS_OF_PARTICLES;
	LaGenMatDouble afull, nfull;
	num=1;//At this stage, I suppose there are only one kind of particle on each site.
	a.resize(num);
	n.resize(num);

	base.read(basefin);
  //std::cout<<base<<std::endl;
	ReadOneOperator(siteopfin, afull);
  //std::cout<<afull<<std::endl;
	ReadOneOperator(siteopfin, nfull);
  //std::cout<<nfull<<std::endl;
	if(base.subnum==1)
		a[0]=Rmatrix(afull,base,base,1);
	else
		a[0]=Rmatrix(afull,base,base,3);
	//std::cout<<a[0]<<std::endl;
	n[0]=Rmatrix(nfull,base,base,1);
	//std::cout<<n[0]<<std::endl;
}


/*!
    \fn Site::write(std::ofstream &fout)
 */
void Site::write(std::ofstream &fout)
{
  fout.write(&value_type,sizeof value_type);
  base.write(fout);
  for(int i=0;i<num;i++)
  {
  if(value_type=='r')
  {
    //std::cout<<a[i]<<std::endl;
    a[i].write(fout);
    n[i].write(fout);
  }
  else
  {
    aC[i].write(fout);
    nC[i].write(fout);
  }
  }
}


/*!
    \fn Site::read(std::ifstream &fin)
 */
void Site::read(std::ifstream &fin)
{
  fin.read(&value_type,sizeof value_type);
  base.read(fin);
  num=NUMBER_OF_KINDS_OF_PARTICLES;
  a.resize(num);
  n.resize(num);
  for(int i=0;i<num;i++)
  {
  if(value_type=='r')
  {
    a[i].read(fin);
    n[i].read(fin);
  }
  else
  {
    aC[i].read(fin);
    nC[i].read(fin);
  }
  }
}


/*!
    \fn Site::Site(std::ifstream &fin)
 */
Site::Site(std::ifstream &fin)
{
  num=NUMBER_OF_KINDS_OF_PARTICLES;
  read(fin);
}


/*!
    \fn Site::operator=(Site &s)
 */

Site& Site::operator=(const Site& s)
{
  value_type=s.value_type;
  num=s.num;
    a=s.a;
    n=s.n;

    aC=s.aC;
    nC=s.nC;

  base=s.base;

}



/*!
    \fn Site::toComplex()
 */
void Site::toComplex()
{
  if(value_type=='r')
  {
    value_type='c';
    aC.resize(num);
    nC.resize(num);
    for(int i=0;i<num;i++)
    {
      R2C(a[i],aC[i]);
      R2C(n[i],nC[i]);
      a[i].delmat();
      n[i].delmat();
    }
  }
}

std::ostream & operator<<(std::ostream& os, const Site& site)
{
  std::cout<<"===Site Information==="<<std::endl;
  std::cout<<site.base<<std::endl;
  for(int i=0;i<site.num;i++)
  {
  std::cout<<"Particle "<<i<<":"<<std::endl;
    std::cout<<site.a[i]<<std::endl;
    std::cout<<site.n[i]<<std::endl;
  }
  return os;
}


/*!
\fn Site::multsignmat()
 */
void Site::multsignmat()
{
  for(int i=0;i<num;i++)
  {
    Rmatrix signmat;
    signmat.gensignmat(base,i);
    // std::cout<<a[i]<<std::endl;
    //std::cout<<signmat<<std::endl;
    a[i]=a[i]*signmat;
  }
}



/*!
    \fn Site::genfermion()
 */
/*
void Site::genfermion()
{
  int dim=4;
  base.Dim=4;
  base.subnum=4;

  base.TargetGQN.resize(dim);
  base.dim.resize(dim);
  base.TargetGQN[0].gqn[0]=0;
  base.TargetGQN[0].gqn[1]=0;
  base.TargetGQN[1].gqn[0]=1;
  base.TargetGQN[1].gqn[1]=-1;
  base.TargetGQN[2].gqn[0]=1;
  base.TargetGQN[2].gqn[1]=1;
  base.TargetGQN[3].gqn[0]=2;
  base.TargetGQN[3].gqn[1]=0;

  base.dim[0]=1;
  base.dim[1]=1;
  base.dim[2]=1;
  base.dim[3]=1;


  a.resize(num);
  n.resize(num);

  for(int i=0;i<num;i++)
  {
    a[i].resize(base,base);
    n[i].resize(base,base);
    a[i].subnum=2;
    a[i].submat.resize(a[i].subnum);
    n[i].subnum=4;///Though there is only two nozero elements, the zeros on the diagnal have their physical meanig.
    n[i].submat.resize(n[i].subnum);
    for(int j=0;j<2;j++)
    {
      a[i].submat[j].resize(1,1);
      n[i].submat[j].resize(1,1);
      n[i].submat[2+j].resize(1,1);
    }
  }

  //a_down
  a[0].pmat(0,1)=0;
  a[0].pmat(2,3)=1;
  //a_up
  a[1].pmat(0,2)=0;
  a[1].pmat(1,3)=1;

  n[0].pmat(0,0)=0;
  n[0].pmat(1,1)=1;
  n[0].pmat(2,2)=2;
  n[0].pmat(3,3)=3;
  n[1].pmat=n[0].pmat;

  a[0].submat[0]=1;
  a[0].submat[1]=1;
  a[1].submat[0]=1;
  a[1].submat[1]=-1;
  n[0].submat[0]=0;
  n[0].submat[1]=1;
  n[0].submat[2]=0;
  n[0].submat[3]=1;
  n[1].submat[0]=0;
  n[1].submat[1]=0;
  n[1].submat[2]=1;
  n[1].submat[3]=1;


}
*/

/*!
    \fn Site::genspin()
 */
/*
void Site::genspin()
{
  int dim=2;
  base.Dim=2;
  base.subnum=2;

  base.TargetGQN.resize(dim);
  base.dim.resize(dim);
  base.TargetGQN[0].gqn[0]=-1;
  base.TargetGQN[1].gqn[0]=1;

  base.dim[0]=1;
  base.dim[1]=1;

  a.resize(num);
  n.resize(num);


  a[0].resize(base,base);
  a[0].pmat(0,1)=0;
  a[0].subnum=1;
  a[0].submat.resize(a[0].subnum);
  a[0].submat[0].resize(1,1);
  a[0].submat[0]=1;

  n[0].resize(base,base);
  n[0].pmat(0,0)=0;
  n[0].pmat(1,1)=1;
  n[0].subnum=2;
  n[0].submat.resize(n[0].subnum);
  n[0].submat[0].resize(1,1);
  n[0].submat[0]=-0.5;
  n[0].submat[1].resize(1,1);
  n[0].submat[1]=0.5;

  /*
#if USE_EVOLVE && !USE_INFINITE_DMRG

#endif
  */
//}


/*!
    \fn Site::genspinlessfermion()
 */
/*
void Site::genspinlessfermion()
{
  int dim=2;
  base.Dim=2;
  base.subnum=2;

  base.TargetGQN.resize(dim);
  base.dim.resize(dim);
  base.TargetGQN[0].gqn[0]=0;
  base.TargetGQN[1].gqn[0]=1;

  base.dim[0]=1;
  base.dim[1]=1;

  a.resize(num);
  n.resize(num);


  a[0].resize(base,base);
  a[0].pmat(0,1)=0;
  a[0].subnum=1;
  a[0].submat.resize(a[0].subnum);
  a[0].submat[0].resize(1,1);
  a[0].submat[0]=1;

  n[0].resize(base,base);
  n[0].pmat(0,0)=0;
  n[0].pmat(1,1)=1;
  n[0].subnum=2;
  n[0].submat.resize(n[0].subnum);
  n[0].submat[0].resize(1,1);
  n[0].submat[0]=0;
  n[0].submat[1].resize(1,1);
  n[0].submat[1]=1;
}

*/
/*!
    \fn Site::genspinFTDMRG()
 */
/*
void Site::genspinFTDMRG()
{
  LaGenMatDouble sa(2,2);
  LaGenMatDouble sc(2,2);
  LaGenMatDouble sz(2,2);
  LaGenMatDouble u(4,4);
  LaGenMatDouble i2;
  LaGenMatDouble afull,cfull,nfull;

  sa=0;
  sc=0;
  sz=0;
  sa(0,1)=1;
  sc(1,0)=1;

  sz(0,0)=0;
  sz(1,1)=1;

  i2=i2.eye(2,2);
  afull=kron(sa,i2);
  cfull=kron(sc,i2);
  nfull=kron(sz,i2);


  u=0;
  u(0,0)=1;
  u(3,3)=1;
  u(1,1)=1/sqrt(2.);
  u(1,2)=1/sqrt(2.);
  u(2,1)=1/sqrt(2.);
  u(2,2)=-1/sqrt(2.);
  u(3,3)=1;

  LaGenMatDouble temp(4,4);
  Blas_Mat_Mat_Trans_Mult(afull,u,temp);
  Blas_Mat_Mat_Mult(u,temp,afull);

  Blas_Mat_Mat_Trans_Mult(cfull,u,temp);
  Blas_Mat_Mat_Mult(u,temp,cfull);

  Blas_Mat_Mat_Trans_Mult(nfull,u,temp);
  Blas_Mat_Mat_Mult(u,temp,nfull);

  base.Dim=4;
  base.subnum=3;
  base.dim.resize(base.subnum);
  base.TargetGQN.resize(base.subnum);
  base.dim[0]=1;
  base.dim[1]=2;
  base.dim[2]=1;
  base.TargetGQN[0].gqn[0]=-1;
  base.TargetGQN[1].gqn[0]=0;
  base.TargetGQN[2].gqn[0]=1;

  a.resize(num);
  c.resize(num);
  n.resize(num);
  GQN gqn;
  a[0]=Rmatrix(afull,base,base,3,gqn);
  //std::cout<<a;
  c[0]=Rmatrix(cfull,base,base,4,gqn);
  //std::cout<<c;
  n[0]=Rmatrix(nfull,base,base,1,gqn);
  //std::cout<<n;

}
*/

