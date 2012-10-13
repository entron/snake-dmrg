#ifndef GQNMAT_H
#define GQNMAT_H

/**
Class of partitioned matrix accoding to good quantum number.

@author Cheng Guo & Hu Shijie
*/

#include "gqnbase.h"
#include "public.h"


#include <gmi.h>
#include <blas3pp.h>
#include <blaspp.h>
#include <algorithm>
#include <vector>
#include <iostream>
using namespace std;


template <class MType>
class GQNMat{
public:
  typedef MType value_type;

  ///Number of sub matrix
  int subnum;

  ///column an row base
  GQNBase colbase,rowbase;

  ///Sub matrix
  vector<MType> submat;

  ///Position matrix. submat[pmat(i,j)] is the (i,j) submatrix
  LaGenMatInt pmat;

  LaGenMatInt temppmat;
  int tempsubnum;
  vector<MType> tempsubmat;

public:
  GQNMat();
  ~GQNMat();
  GQNMat(ifstream& fin);
  GQNMat(const GQNMat& mat);
  GQNMat(const GQNBase &row,const GQNBase &col);

  GQNMat(MType& mat, GQNBase &rowbase,GQNBase &colbase,int mode);
  GQNMat(MType& mat, GQNBase &rowbase,GQNBase &colbase,vector<GQN> &TargetGQN);
  void inject(const GQNMat& mat);
  //  ostream& operator<<(ostream& os, const GQNMat& mat);
  void write(ofstream& fout);
  void read(ifstream& fin);

  void scale(typename MType::value_type s);

  ///Generate eye matrix with the same base as mat
  void geneye(const GQNMat& mat);
  void geneye(const GQNBase &base);

  ///Generate sign matrix for fermion
  void gensignmat(const GQNBase &base,int type);

  void resize(const GQNBase &row,const GQNBase &col);

  void normalize();
  void delmat();

  typename MType::value_type trace();

  void Convert2Full(MType & fullmat);

  GQNMat& operator=(const GQNMat&);
  void operator+=(const GQNMat& mat);
  GQNMat<MType>& operator=(double s);
  GQNMat<MType> operator*(const GQNMat<MType>&) const;
  GQNMat<MType> operator+(const GQNMat<MType> &b) const;

  void order();
  void combine();

private:
    ///Determin whether the storage postion is equal.
  int pmatequal(const GQNMat& mat) const;

};



template<class MType>
GQNMat<MType>::GQNMat()
{
  subnum=0;
  tempsubnum=0;
}

template<class MType>
GQNMat<MType>::~GQNMat()
{
}


/*!
\fn GQNMat::GQNMat(ifstream& fin)
 */
template<class MType>
GQNMat<MType>::GQNMat(ifstream& fin)
{
  read(fin);
}


/*!
\fn GQNMat::GQNMat(GQNMat& mat)
 */
template<class MType>
GQNMat<MType>::GQNMat(const GQNMat& mat)
{
  *this=mat;
}

template<class MType>
GQNMat<MType>::GQNMat(const GQNBase &row,const GQNBase &col)
{
  rowbase=row;
  colbase=col;
  pmat.resize(rowbase.subnum,colbase.subnum);
  pmat=-1;
  subnum=0;
  tempsubnum=0;
}

/*!
\fn GQNMat::GQNMat(MType& mat, GQNBase &rowbase,GQNBase &colbase,int mode)
 */
///Convert full matrix to partitioned matrix
template<class MType>
GQNMat<MType>::GQNMat(MType& mat, GQNBase &row,GQNBase &col,int mode)
{
  subnum=0;
  tempsubnum=0;
  rowbase=row;
  colbase=col;
  pmat.resize(rowbase.subnum,colbase.subnum);
  pmat=-1;
  int rowstart, colstart;
  GQN tempTargetGQN;

  switch(mode)
  {
  case 1:///Diagnal
    colstart=0;
    for(int j=0;j<colbase.subnum;j++)
    {
      rowstart=0;
      for(int i=0;i<rowbase.subnum;i++)
      {
        if(rowbase.subgqn[i]==colbase.subgqn[j])
        {
          subnum++;
          submat.resize(subnum);
          pmat(i,j)=subnum-1;
          submat[subnum-1]=mat(LaIndex(rowstart,rowstart+rowbase.dim[i]-1),LaIndex(colstart,colstart+colbase.dim[j]-1));
          break;
        }
        rowstart+=rowbase.dim[i];
      }
      colstart+=colbase.dim[j];
    }
    break;
  case 2:///Antidiagnal
    colstart=0;
    tempTargetGQN=colbase.subgqn[0]+rowbase.subgqn[rowbase.subnum-1];
    for(int j=0;j<colbase.subnum;j++)
    {
      rowstart=0;
      for(int i=0;i<rowbase.subnum;i++)
      {
        if(rowbase.subgqn[i]+colbase.subgqn[j]==tempTargetGQN)
        {
          subnum++;
          submat.resize(subnum);
          pmat(i,j)=subnum-1;
          submat[subnum-1]=mat(LaIndex(rowstart,rowstart+rowbase.dim[i]-1),LaIndex(colstart,colstart+colbase.dim[j]-1));
          break;
        }
        rowstart+=rowbase.dim[i];
      }
      colstart+=colbase.dim[j];
    }
    break;
  case 3:///Upper diagnal
    colstart=0;
    tempTargetGQN=colbase.subgqn[1]-rowbase.subgqn[0];
    for(int j=0;j<colbase.subnum;j++)
    {
      rowstart=0;
      for(int i=0;i<rowbase.subnum;i++)
      {
        if(rowbase.subgqn[i]+tempTargetGQN==colbase.subgqn[j])
        {
          subnum++;
          submat.resize(subnum);
          pmat(i,j)=subnum-1;


          cout<<rowstart+rowbase.dim[i]-1<<endl;
          cout<<colstart+colbase.dim[j]-1<<endl;
          submat[subnum-1]=mat(LaIndex(rowstart,rowstart+rowbase.dim[i]-1),LaIndex(colstart,colstart+colbase.dim[j]-1));
          break;
        }
        //cout<<rowstart<<endl;
        rowstart+=rowbase.dim[i];
        //cout<<rowstart<<endl;
      }
      //cout<<colstart<<endl;
      colstart+=colbase.dim[j];
      //cout<<colstart<<endl;
    }
    break;
  case 4:///lower diagnal
    colstart=0;
    tempTargetGQN=rowbase.subgqn[1]-colbase.subgqn[0];
    for(int j=0;j<colbase.subnum;j++)
    {
      rowstart=0;
      for(int i=0;i<rowbase.subnum;i++)
      {
        if(rowbase.subgqn[i]-tempTargetGQN==colbase.subgqn[j])
        {
          subnum++;
          submat.resize(subnum);
          pmat(i,j)=subnum-1;
          submat[subnum-1]=mat(LaIndex(rowstart,rowstart+rowbase.dim[i]-1),LaIndex(colstart,colstart+colbase.dim[j]-1));
          break;
        }
        rowstart+=rowbase.dim[i];
      }
      colstart+=colbase.dim[j];
    }
    break;
  }
}


template<class MType>
GQNMat<MType>::GQNMat(MType& mat, GQNBase &row,GQNBase &col,vector<GQN> &TargetGQN)
{
  subnum=0;
  tempsubnum=0;
  rowbase=row;
  colbase=col;
  pmat.resize(rowbase.subnum,colbase.subnum);
  pmat=-1;
  int rowstart,colstart;


    colstart=0;
    for(int j=0;j<colbase.subnum;j++)
    {
      rowstart=0;
      for(int i=0;i<rowbase.subnum;i++)
      {
        if(rowbase.subgqn[i]+colbase.subgqn[j]==TargetGQN)
        {
          subnum++;
          submat.resize(subnum);
          pmat(i,j)=subnum-1;
          submat[subnum-1]=mat(LaIndex(rowstart,rowstart+rowbase.dim[i]-1),LaIndex(colstart,colstart+colbase.dim[j]-1));
          //break;
        }
        rowstart+=rowbase.dim[i];
      }
      colstart+=colbase.dim[j];
    }


}


/*!
\fn GQNMat::inject(GQNMat& mat)
 */
template<class MType>
void GQNMat<MType>::inject(const GQNMat<MType>& mat)
{
  if(pmatequal(mat))
    for(int i=0;i<rowbase.subnum;i++)
      for(int j=0;j<colbase.subnum;j++)
        if(pmat(i,j)>=0)
          submat[pmat(i,j)].inject(mat.submat[pmat(i,j)]);
}


/*!
\fn GQNMat::write(ofstream& fout)
 */
template<class MType>
void GQNMat<MType>::write(ofstream& fout)
{
  fout.write((char*)&subnum,sizeof subnum);
  rowbase.write(fout);
  colbase.write(fout);
  writemat(fout,pmat);
  for(int i=0;i<subnum;i++)
    writemat(fout,submat[i]);
}


/*!
\fn GQNMat::read(ifstream& fin)
 */
template<class MType>
void GQNMat<MType>::read(ifstream& fin)
{
  fin.read((char*)&subnum,sizeof subnum);
  rowbase.read(fin);
  colbase.read(fin);
  readmat(fin,pmat);
  submat.resize(subnum);
  for(int i=0;i<subnum;i++)
    readmat(fin,submat[i]);
}



/*!
\fn GQNMat::scale(MType::value_type s)
 */
template<class MType>
void GQNMat<MType>::scale(typename MType::value_type s)
{
  for(int i=0;i<subnum;i++)
    submat[i].scale(s);
}


/*!
\fn GQNMat::pmatequal(GQNMat& mat)
 */
template<class MType>
int GQNMat<MType>::pmatequal(const GQNMat<MType>& mat) const
{
  if(colbase!=mat.colbase||rowbase!=mat.rowbase||subnum!=mat.subnum)
  {
    cout<<"Unequal base or subnum!"<<endl;
    return 0;
  }

  for(int i=0;i<rowbase.subnum;i++)
    for(int j=0;j<colbase.subnum;j++)
      if(pmat(i,j)*mat.pmat(i,j)<0)
        return 0;

  return 1;
}


/*!
\fn GQNMat::geneye(GQNMat& mat)
 */
template<class MType>
void GQNMat<MType>::geneye(const GQNMat<MType>& mat)
{
  if(mat.rowbase!=mat.colbase)
  {
    cout<<"Different row and column base. Cannot generate eye matrix!"<<endl;
    return;
  }
  else
  {
    subnum=mat.rowbase.subnum;
    submat.resize(subnum);
    rowbase=mat.rowbase;
    colbase=rowbase;
    pmat.resize(rowbase.subnum,rowbase.subnum);
    pmat=-1;
    for(int i=0;i<rowbase.subnum;i++)
    {
      pmat(i,i)=i;
      submat[i]=submat[i].eye(rowbase.dim[i]);
    }
  }
}

template<class MType>
void GQNMat<MType>::geneye(const GQNBase &base)
{
  rowbase=base;
  colbase=base;
  subnum=rowbase.subnum;
  submat.resize(subnum);
  pmat.resize(rowbase.subnum,rowbase.subnum);
  pmat=-1;
  for(int i=0;i<rowbase.subnum;i++)
  {
    pmat(i,i)=i;
    submat[i]=submat[i].eye(rowbase.dim[i]);
  }
}


template<class MType>
void GQNMat<MType>::gensignmat(const GQNBase &base,int type)
{
  rowbase=base;
  colbase=base;
  subnum=rowbase.subnum;
  submat.resize(subnum);
  pmat.resize(rowbase.subnum,rowbase.subnum);
  pmat=-1;
  if(type==0)//down spin
    for(int i=0;i<rowbase.subnum;i++)
    {
      pmat(i,i)=i;
      submat[i].resize(1,1);
      if(i<2)    submat[i]=1;
      else submat[i]=-1;
    }
  else//up spin
    for(int i=0;i<rowbase.subnum;i++)
    {
      pmat(i,i)=i;
      submat[i].resize(1,1);
      if(i==0||i==2)
        submat[i]=1;
      else
        submat[i]=-1;
    }
}


template<class MType>
void GQNMat<MType>::Convert2Full(MType & fullmat)
{
  fullmat.resize(rowbase.Dim,colbase.Dim);
  fullmat=0;
  int rowstart,colstart;
  rowstart=0;
  for(int i=0;i<rowbase.subnum;i++)
  {
    colstart=0;
    for(int j=0;j<colbase.subnum;j++)
    {
      if(pmat(i,j)!=-1)
        fullmat(LaIndex(rowstart,rowstart+rowbase.dim[i]-1),LaIndex(colstart,colstart+colbase.dim[j]-1)).inject(submat[pmat(i,j)]);
      colstart+=colbase.dim[j];
    }
    rowstart+=rowbase.dim[i];
  }
}


template<class MType>
void GQNMat<MType>::resize(const GQNBase &row,const GQNBase &col)
{
  rowbase=row;
  colbase=col;
  pmat.resize(rowbase.subnum,colbase.subnum);
  pmat=-1;
  subnum=0;
  tempsubnum=0;
}

template<class MType>
typename MType::value_type GQNMat<MType>::trace()
{
  ///This function is do not write in a easy to understand way because of lapackpp don't provied me with easy to use functions and operators to realize is function.Say there are no operators for COMPLEX add.
  int isfirst=0;
  MType t;
  t.resize(1,1);
  if(rowbase!=colbase)
  {
    cout<<"rowbase!=colbase"<<endl;
  }
  else
  {
    for(int i=0;i<rowbase.subnum;i++)
      if(pmat(i,i)!=-1)
        if(isfirst==0)
        {
          t=submat[pmat(i,i)].trace();
          isfirst=1;
        }
        else
          t+=submat[pmat(i,i)].trace();
        return t(0,0);
  }
}


template<class MType>
void GQNMat<MType>::normalize()
{
  double normsq=0,subnorm,norm;
  for(int i=0;i<rowbase.subnum;i++)
    for(int j=0;j<colbase.subnum;j++)
      if(pmat(i,j)!=-1)
      {
        subnorm=Blas_NormF(submat[pmat(i,j)]);
        normsq+=subnorm*subnorm;
      }
  norm=sqrt(normsq);
  for(int i=0;i<subnum;i++)
    submat[i].scale(1/norm);
}


template<class MType>
void GQNMat<MType>::delmat()
{
  subnum=0;
  tempsubnum=0;
  pmat.resize(0,0);
  temppmat.resize(0,0);
  submat.resize(0);
  tempsubmat.resize(0);

}


///IMPORTANT!!! c shoulb be an empoty matrix that is pmat=-1 subnum=0 before using these
///multiplication functions
template<class MType>
void Mat_Mat_Mult(const GQNMat<MType> &a,const GQNMat<MType> &b,GQNMat<MType> &c)
{
  if(a.colbase==b.rowbase)
  {
    if(a.rowbase==c.rowbase && b.colbase==c.colbase)
    {
      for(int j=0;j<c.colbase.subnum;j++)
        for(int i=0;i<c.rowbase.subnum;i++)
          for(int k=0;k<a.colbase.subnum;k++)
            if(a.pmat(i,k)>=0 && b.pmat(k,j)>=0)

              if(c.pmat(i,j)>=0)
                Blas_Mat_Mat_Mult(a.submat[a.pmat(i,k)],b.submat[b.pmat(k,j)],c.submat[c.pmat(i,j)],1.0,1.0);
      else
      {
        c.subnum++;
        c.submat.resize(c.subnum);
        c.pmat(i,j)=c.subnum-1;
        c.submat[c.subnum-1].resize(a.submat[a.pmat(i,k)].size(0),b.submat[b.pmat(k,j)].size(1));
        Blas_Mat_Mat_Mult(a.submat[a.pmat(i,k)],b.submat[b.pmat(k,j)],c.submat[c.pmat(i,j)]);
      }
    }
    else cout<<"a.rowbase!=c.rowbase || b.colbase!=c.colbase"<<endl;
  }
  else cout<<"a.colbase!=b.rowbase"<<endl;
}


template<class MType>
void Mat_Trans_Mat_Mult(const GQNMat<MType> &a,const GQNMat<MType> &b,GQNMat<MType> &c)
{
  if(a.rowbase==b.rowbase)
  {
    if(a.colbase==c.rowbase && b.colbase==c.colbase)
    {
      for(int j=0;j<c.colbase.subnum;j++)
        for(int i=0;i<c.rowbase.subnum;i++)
          for(int k=0;k<a.rowbase.subnum;k++)
            if(a.pmat(k,i)>=0 && b.pmat(k,j)>=0)

              if(c.pmat(i,j)>=0)
                Blas_Mat_Trans_Mat_Mult(a.submat[a.pmat(k,i)],b.submat[b.pmat(k,j)],c.submat[c.pmat(i,j)],1.0,1.0);
      else
      {
        c.subnum++;
        c.submat.resize(c.subnum);
        c.pmat(i,j)=c.subnum-1;
        c.submat[c.subnum-1].resize(a.submat[a.pmat(k,i)].size(1),b.submat[b.pmat(k,j)].size(1));
        Blas_Mat_Trans_Mat_Mult(a.submat[a.pmat(k,i)],b.submat[b.pmat(k,j)],c.submat[c.pmat(i,j)]);
        //cout<<c.submat[c.pmat(i,j)]<<endl;
      }
    }
    else cout<<"a.colbase!=c.rowbase || b.colbase!=c.colbase"<<endl;
  }
  else cout<<"a.rowbase!=b.rowbase"<<endl;
}


template<class MType>
void Mat_Mat_Trans_Mult(const GQNMat<MType> &a,const GQNMat<MType> &b,GQNMat<MType> &c)
{
  if(a.colbase==b.colbase)
  {
    if(a.rowbase==c.rowbase && b.rowbase==c.colbase)
    {
      for(int j=0;j<c.colbase.subnum;j++)
        for(int i=0;i<c.rowbase.subnum;i++)
          for(int k=0;k<a.colbase.subnum;k++)
            if(a.pmat(i,k)>=0 && b.pmat(j,k)>=0)

              if(c.pmat(i,j)>=0)
                Blas_Mat_Mat_Trans_Mult(a.submat[a.pmat(i,k)],b.submat[b.pmat(j,k)],c.submat[c.pmat(i,j)],1.0,1.0);
      else
      {
        c.subnum++;
        c.submat.resize(c.subnum);
        c.pmat(i,j)=c.subnum-1;
        c.submat[c.subnum-1].resize(a.submat[a.pmat(i,k)].size(0),b.submat[b.pmat(j,k)].size(0));
        Blas_Mat_Mat_Trans_Mult(a.submat[a.pmat(i,k)],b.submat[b.pmat(j,k)],c.submat[c.pmat(i,j)]);
      }
    }
    else cout<<"a.rowbase!=c.rowbase || b.rowbase!=c.colbase"<<endl;
  }
  else cout<<"a.colbase!=b.colbase"<<endl;
}

template<class MType>
void Transpose(const GQNMat<MType> &a,GQNMat<MType> &b)
{
  //Base
  b.colbase=a.rowbase;
  b.rowbase=a.colbase;
  b.subnum=a.subnum;
  b.submat.resize(b.subnum);
  b.pmat.resize(b.rowbase.subnum,b.colbase.subnum);
  //pmat
  for(int i=0;i<a.rowbase.subnum;i++)
    for(int j=0;j<a.colbase.subnum;j++)
      b.pmat(j,i)=a.pmat(i,j);
  //submat
  for(int i=0;i<b.subnum;i++)
  {
    b.submat[i].resize(a.submat[i].size(1),a.submat[i].size(0));
    for(int r=0;r<a.submat[i].size(0);r++)
      for(int c=0;c<a.submat[i].size(1);c++)
        b.submat[i](c,r)=a.submat[i](r,c);
    /*
    MType eyemat;
    cout<<a.submat[i]<<endl;
    cout<<a.submat[i].size(0)<<endl;
    cout<<eyemat.eye(a.submat[i].size(0),a.submat[i].size(0))<<endl;
    eyemat=eyemat.eye(a.submat[i].size(0),a.submat[i].size(0));
    b.submat[i].resize(a.submat[i].size(1),a.submat[i].size(0));
    //cout<<eyemat<<endl;
    //cout<<b.submat[i]<<endl;
    Blas_Mat_Trans_Mat_Mult(a.submat[i],eyemat,b.submat[i]);
    */
  }
}



///Caculate two matrix dot product.a_{i,j}*b_{i,j}.
template<class MType>
typename MType::value_type Mat_Dot_Prod(const GQNMat<MType> &a,const GQNMat<MType> &b)
{
  typename MType::value_type prod=0;
  if(a.rowbase==b.rowbase && a.colbase==b.colbase)
  {
    for(int i=0;i<a.rowbase.subnum;i++)
      for(int j=0;j<a.colbase.subnum;j++)
        if(a.pmat(i,j)!=-1 && b.pmat(i,j)!=-1)
        {
          for(int m=0;m<a.submat[a.pmat(i,j)].size(0);m++)
            for(int n=0;n<a.submat[a.pmat(i,j)].size(1);n++)
                prod+=a.submat[a.pmat(i,j)](m,n)*b.submat[b.pmat(i,j)](m,n);
        }
  }
  else
    cout<<"Matrix bases are not equal, can't caculate Mat_Dot_Prod"<<endl;
  return prod;

}

template<class MType>
void kron(GQNMat<MType> &a,GQNMat<MType> &b,GQNMat<MType> &c)
{
  int arow=a.rowbase.subnum;
  int acol=a.colbase.subnum;
  int brow=b.rowbase.subnum;
  int bcol=b.colbase.subnum;
  int p,q;

  c.rowbase=kron(a.rowbase,b.rowbase);
  // cout<<a.rowbase<<endl;
  //cout<<b.rowbase<<endl;
  //cout<<c.rowbase<<endl;
  c.colbase=kron(a.rowbase,b.rowbase);
  c.temppmat.resize(arow*brow,acol*bcol);
  c.temppmat=-1;
  c.tempsubnum=0;
  ///Direct product(http://mathworld.wolfram.com/MatrixDirectProduct.html)
  ///c_{pq}=a_{ij}*b_{kl}
  ///p=brow*(i-1)+k
  ///q=bcol*(j-1)+l
  for(int i=0;i<arow;i++)
    for(int j=0;j<acol;j++)
      for(int k=0;k<brow;k++)
        for(int l=0;l<bcol;l++)
        {
          p=brow*i+k;
          q=bcol*j+l;
          if(a.pmat(i,j)>=0 && b.pmat(k,l)>=0)
          {
            c.tempsubnum++;
            c.tempsubmat.resize(c.tempsubnum);
            c.temppmat(p,q)=c.tempsubnum-1;
            c.tempsubmat[c.tempsubnum-1]=kron(a.submat[a.pmat(i,j)],b.submat[b.pmat(k,l)]);
          }
        }
  c.order();
  c.combine();

  c.temppmat.resize(0,0);
  c.tempsubmat.resize(0);
  c.tempsubnum=0;

}

///Change pmat so that the matrix base is aligned according to good quantum number.
template<class MType>
void GQNMat<MType>::order()
{
  ///Order pmat
  //cout<<temppmat<<endl;
  ordermat(rowbase.map,colbase.map,temppmat);
  //cout<<temppmat<<endl;
  ///Order tempdim
  vector<int> v;
  v=rowbase.tempdim;
  for(int i=0;i<v.size();i++)
    rowbase.tempdim[i]=v[rowbase.map[i]];
  v=colbase.tempdim;
  for(int i=0;i<v.size();i++)
    colbase.tempdim[i]=v[colbase.map[i]];
}


///Combine the submatrix which belong to the same good quantum number
template<class MType>
void GQNMat<MType>::combine()
{
  //cout<<rowbase.tempTargetGQN.size()<<endl;
  vector<GQN> rowsubgqn(rowbase.tempsubgqn);
  vector<GQN> colsubgqn(colbase.tempsubgqn);
  pmat.resize(rowbase.subnum,colbase.subnum);
  pmat=-1;

  sort(rowsubgqn.begin(),rowsubgqn.end());
  sort(colsubgqn.begin(),colsubgqn.end());

  int rowstart,colstart;
  for(int n=0;n<colbase.subnum;n++)
    for(int m=0;m<rowbase.subnum;m++)
    {

      colstart=0;
      for(int j=0;j<colsubgqn.size();j++)
      {
        if(colsubgqn[j]!=colbase.subgqn[n]) continue;
        else
        {
          rowstart=0;
          for(int i=0;i<rowsubgqn.size();i++)
          {
            if(rowsubgqn[i]!=rowbase.subgqn[m]) continue;
            else
            {
              if(temppmat(i,j)==-1)
              {
                rowstart+=rowbase.tempdim[i];
                continue;
              }
              if(pmat(m,n)==-1)
              {
                subnum++;
                submat.resize(subnum);
                pmat(m,n)=subnum-1;
                submat[subnum-1].resize(rowbase.dim[m],colbase.dim[n]);
                submat[subnum-1]=0;
              }
              submat[subnum-1](LaIndex(rowstart,rowstart+rowbase.tempdim[i]-1),LaIndex(colstart,colstart+colbase.tempdim[j]-1)).inject(tempsubmat[temppmat(i,j)]);
              rowstart+=rowbase.tempdim[i];
            }
          }
          colstart+=colbase.tempdim[j];
        }
      }
    }
}

template<class MType>
ostream & operator<<(ostream& os, const GQNMat<MType>& mat)
{
  cout<<"+++++GQNMat Information+++++"<<endl;
  cout<<"Subnum="<<mat.subnum<<endl;
cout<<"Postion matrix is: "<<endl;
  cout<<mat.pmat<<endl;
  for(int i=0;i<mat.subnum;i++)
  {
  cout<<"The "<<i<<"th submatrix is: "<<endl;
    cout<<mat.submat[i]<<endl;
  }
  cout<<"+++++End of GQNMat Information+++++"<<endl;
  return os;

}

/*!
\fn GQNMat::operator=(GQNMat& mat)
 */

template<class MType>
GQNMat<MType>& GQNMat<MType>::operator=(const GQNMat<MType> &mat)
{
  subnum=mat.subnum;
  colbase=mat.colbase;
  rowbase=mat.rowbase;
  pmat=mat.pmat;
  submat=mat.submat;
  return *this;
}

/*!
\fn GQNMat::operator+=(GQNMat& mat)
 */
template<class MType>
void GQNMat<MType>::operator+=(const GQNMat<MType>& mat)
{
  if(rowbase==mat.rowbase&&colbase==mat.colbase)
  {
    for(int i=0;i<rowbase.subnum;i++)
      for(int j=0;j<colbase.subnum;j++)
      {
        if(pmat(i,j)>= 0&&mat.pmat(i,j)>= 0)
        {
           addmat(submat[pmat(i,j)],mat.submat[mat.pmat(i,j)]);
          //submat[pmat(i,j)]=submat[pmat(i,j)]+mat.submat[mat.pmat(i,j)];///Need change for efficiency
        }
        else if(pmat(i,j)==-1&&mat.pmat(i,j)>= 0)
        {
          subnum++;
          submat.resize(subnum);
          pmat(i,j)=subnum-1;
          submat[subnum-1]=mat.submat[mat.pmat(i,j)];
        }
      }
  }
  else
  {
    cout<<"Matrix size are unequal. Can not add!"<<endl;
  }
}


/*!
\fn GQNMat::operator=(MType::value_type s)
 */

template<class MType>
GQNMat<MType>& GQNMat<MType>::operator=(double s)
{
  for(int i=0;i<rowbase.subnum;i++)
    for(int j=0;j<colbase.subnum;j++)
    {
      if(pmat(i,j)<0) continue;
      if(submat[pmat(i,j)].size(0)!=rowbase.dim[i]||submat[pmat(i,j)].size(1)!=colbase.dim[j])
      {
        submat[pmat(i,j)].resize(rowbase.dim[i],colbase.dim[j]);
      }
      submat[pmat(i,j)]=s;
    }
}


template<class MType>
GQNMat<MType> GQNMat<MType>::operator *(const GQNMat<MType> &b) const
{
  GQNMat<MType> newmat(rowbase,b.colbase);
  Mat_Mat_Mult(*this,b,newmat);
  return newmat;
}

template<class MType>
GQNMat<MType> GQNMat<MType>::operator+(const GQNMat<MType> &b) const
{
  GQNMat<MType> r;
  r=*this;
  r+=b;
  return r;
}



void R2C(GQNMat<LaGenMatDouble> &r,GQNMat<LaGenMatComplex> &c);
void C2R(GQNMat<LaGenMatComplex> &c,GQNMat<LaGenMatDouble> &r);
#endif
