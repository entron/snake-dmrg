#include "blocham.h"
#include "dtmat.h"
#include "site.h"
#include "public.h"

#include<blas3pp.h>
#include<blaspp.h>

#include <iostream>
using namespace std;

BlocHam::BlocHam()
{
  value_type='r';
}

BlocHam::~BlocHam()
{
}

/*!
\fn BlocHam::BlocHam(Site &first)
 */
BlocHam::BlocHam(Site &first,double OnSiteE)
{
	value_type=first.value_type;
	base=first.base;
	if(OnSiteE!=0)
	{
		cout<<"!!!!!!!!!!!!!!!!!!!!!Speical Case: H=OnSiteE*sigma_x!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		H=first.n[0];
		H.scale(OnSiteE);
  	//cout<<H<<endl;
	}
	else{
		if(value_type=='r')
		{
			H.resize(base,base);
			H=0;
		}
		else
		{
			HC.resize(base,base);
			HC=0;
		}
	}
}

/*!
\fn BlocHam::BlocHam(ifstream &fin)
 */
BlocHam::BlocHam(ifstream &fin)
{
  read(fin);
}




/*!
\fn BlocHam::newblocham(Site &add)
 */

BlocHam::BlocHam(BlocHam *oldham,Site &add,Rmatrix &Hinter,char hand)
{
  value_type='r';
  Rmatrix identity;
  identity.geneye(add.base);
///If the site is added to the right of the block
  if(hand=='r')
  {
    ///Caculate new H
    kron(oldham->H,identity,H);
    H+=Hinter;
    base=kron(oldham->base,add.base);
    //cout<<oldham->base<<endl;
    //cout<<add.base<<endl;
    //cout<<base<<endl;
    // cout<<H<<endl;
  }
  else
  {
    kron(identity,oldham->H,H);
    H+=Hinter;
    base=kron(add.base,oldham->base);
  }
}

BlocHam::BlocHam(BlocHam *oldham,Site &add,Cmatrix &Hinter,char hand)
{
  value_type='c';
  Cmatrix identity;
  identity.geneye(add.base);
  if(hand=='r')
  {
    ///Caculate new H
    kron(oldham->HC,identity,HC);
    HC+=Hinter;
    base=kron(oldham->base,add.base);
  }
  else
  {
    kron(identity,oldham->HC,HC);
    HC+=Hinter;
    base=kron(add.base,oldham->base);
  }
}


/*!
\fn BlocHam::BlocHam(BlocHam *oldham,Site &add,LaGenMatDouble &Hinter,char hand,int NOinter)
 */
BlocHam::BlocHam(BlocHam *oldham,Site &add,char hand)
{
  value_type=oldham->value_type;
  if(value_type=='r')
  {
    Rmatrix identity;
    identity.geneye(add.base);
///If the site is added to the right of the block
    if(hand=='r')
    {
    ///Caculate new H
      kron(oldham->H,identity,H);
      base=kron(oldham->base,add.base);
    }
    else
    {
      kron(identity,oldham->H,H);
      base=kron(add.base,oldham->base);
    }
  }
  else
  {
    Cmatrix identity;
    identity.geneye(add.base);
    if(hand=='r')
    {
    ///Caculate new H
      kron(oldham->HC,identity,HC);
      base=kron(oldham->base,add.base);
    }
    else
    {
      kron(identity,oldham->HC,HC);
      base=kron(add.base,oldham->base);
    }
  }

}


/*!
\fn BlocHam::renorm(LaGenMatDouble &tmat)
 */
void BlocHam::renorm(DTMat &mat)
{
  if(value_type=='r')
  {
    Rmatrix temp(mat.tmatbase,base);
  //cout<<"TranMat"<<endl<<mat.trunmattrans<<endl;
    Mat_Trans_Mat_Mult(mat.trunmat,H,temp);
    H.resize(mat.tmatbase,mat.tmatbase);
    Mat_Mat_Mult(temp,mat.trunmat,H);
  //cout<<"H after renorm is"<<endl<<H<<endl;
  }
  else
  {
    Cmatrix temp(mat.tmatbase,base);
  //cout<<"TranMat"<<endl<<mat.trunmattrans<<endl;
    Mat_Trans_Mat_Mult(mat.trunmatC,HC,temp);
    HC.resize(mat.tmatbase,mat.tmatbase);
    Mat_Mat_Mult(temp,mat.trunmatC,HC);
  //cout<<"H after renorm is"<<endl<<H<<endl;
  }
  base=mat.tmatbase;
}



/*!
    \fn BlocHam::write(ofstream &fout)
 */
void BlocHam::write(ofstream &fout)
{
  fout.write(&value_type,sizeof value_type);
  base.write(fout);
  //for(int i=0;i<base.subnum;i++) writemat(fout,h[i]);
  if(value_type=='r')
  {
    H.write(fout);
  }
  else
  {
    HC.write(fout);
  }

}


/*!
    \fn BlocHam::read(ifstream &fin)
 */
void BlocHam::read(ifstream &fin)
{
  fin.read(&value_type,sizeof value_type);
  base.read(fin);
  //for(int i=0;i<base.subnum;i++) readmat(fin,h[i]);
  //cout<<dim[1]<<endl;
  if(value_type=='r')
  {
    H.read(fin);
  }
  else
  {
    HC.read(fin);
  }
}



/*!
    \fn BlocHam::toComplex()
 */
void BlocHam::toComplex()
{
  value_type='c';
  R2C(H,HC);
  H.delmat();
}