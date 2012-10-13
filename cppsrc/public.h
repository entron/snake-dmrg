#ifndef OPERATION_H
#define OPERATION_H

#include<vector>
#include<fstream>
using namespace std;

#include "gqnbase.h"


#include<gmd.h>
#include<lavd.h>
#define LA_COMPLEX_SUPPORT
#include <gmc.h>
#include <lavc.h>


LaGenMatDouble kron(LaGenMatDouble &a,LaGenMatDouble &b);

LaGenMatComplex kron(LaGenMatComplex &a,LaGenMatComplex &b);

LaGenMatDouble directsum(LaGenMatDouble &a,LaGenMatDouble &b);

void chop(LaGenMatDouble &a,double err=1e-15);

void chop(LaGenMatComplex &a,double err=1e-15);

/**Order matrix bases coording to map
*/
template<class MType>
void ordermat(vector<int> &rowmap,vector<int> &colmap,MType &matrix)
{
  //cout<<matrix<<endl;
  int row=matrix.size(0);
  int col=matrix.size(1);

  MType temp(row,col),tempvec;
  //Order column
  if(colmap.size()!=0) //Do not order column if the colmap is empty.
  {
    for(int i=0;i<col;i++)
    {
      tempvec=matrix(LaIndex(0,row-1),LaIndex(colmap[i]));
      temp(LaIndex(0,row-1),LaIndex(i)).inject(tempvec);
    }
    matrix.inject(temp);
  }
  //cout<<temp<<endl;
  //Order row
  if(rowmap.size()!=0)
  {

    for(int i=0;i<row;i++)
    {
      tempvec=matrix(LaIndex(rowmap[i]),LaIndex(0,col-1));
      temp(LaIndex(i),LaIndex(0,col-1)).inject(tempvec);
    }
    matrix.inject(temp);
  }
  //  cout<<matrix<<endl;
};


template<class MType>
void unordermat(vector<int> &rowmap,vector<int> &colmap,MType &matrix)
{
  //cout<<rowmap<<endl;
  // cout<<colmap<<endl;
  int row=matrix.size(0);
  int col=matrix.size(1);
  MType temp(row,col),tempvec;
  //Order column
  if(colmap.size()!=0)
  {
    for(int i=0;i<col;i++)
    {
      tempvec=matrix(LaIndex(0,row-1),LaIndex(i));
      temp(LaIndex(0,row-1),LaIndex(colmap[i])).inject(tempvec);
    }
    matrix.inject(temp);
  }
  //cout<<temp<<endl;
  //Order row
  if(rowmap.size()!=0)
  {
    for(int i=0;i<row;i++)
    {
      tempvec=matrix(LaIndex(i),LaIndex(0,col-1));
      temp(LaIndex(rowmap[i]),LaIndex(0,col-1)).inject(tempvec);
    }
    matrix.inject(temp);
  }
  //cout<<matrix<<endl;
};

/**Read matrix from file
*/
template<class MType>
void readmat(ifstream &fin,MType &mat)
{
  int row,col;
  fin.read((char*)&row,sizeof row);
  fin.read((char*)&col,sizeof col);
  mat.resize(row,col);
  fin.read((char*)mat.addr(),row*col*sizeof(typename MType::value_type));

};

/**Write matrix information to file
*/
template<class MType>
void writemat(ofstream &fout,MType &mat)
{
  int row,col;
  row=mat.size(0);
  col=mat.size(1);
  fout.write((char*)&row,sizeof row);
  fout.write((char*)&col,sizeof col);
  fout.write((char*)mat.addr(),row*col*sizeof(typename MType::value_type));
};

/**Read matrix from file
*/
template<class VType>
void readvec(ifstream &fin,VType &vec)
{
  int row;
  fin.read((char*)&row,sizeof row);
  vec.resize(row,1);
  fin.read((char*)vec.addr(),row*sizeof(typename VType::value_type));

};

/**Write matrix information to file
*/
template<class VType>
void writevec(ofstream &fout,VType &vec)
{
  int row;
  row=vec.size();
  fout.write((char*)&row,sizeof row);
  fout.write((char*)vec.addr(),row*sizeof(typename VType::value_type));
};

template<class Type>
void readvector(ifstream &fin,vector<Type> &v)
{
  int dim;
  dim=v.size();
  if(dim==0)
  {
    fin.read((char*)&dim,sizeof dim);
    v.resize(dim);
  }
  for(int i=0;i<dim;i++)
  {
    fin.read((char*)&v[i],sizeof(Type));
    //cout<<v[i]<<endl;
  }
};

template<class Type>
void writevector(ofstream &fout,vector<Type> &v)
{
  int dim=v.size();
  fout.write((char*)&dim,sizeof dim);
  for(int i=0;i<dim;i++)
    fout.write((char*)&v[i],sizeof(Type));
};



//void combine(char *s1, int ks,char *s2);

double trace(LaGenMatDouble &mat);

template<class MType,class VType>
void mat2vec(VType &v,MType &mat)
{
  int row=mat.size(0);
  int col=mat.size(1);
  if(v.size()!=row*col) v.resize(row*col,1);
  int pointer=0;
  for(int i=0;i<col;i++)
  {
    v(LaIndex(pointer,pointer+row-1)).inject(mat(LaIndex(),LaIndex(i)));
    pointer+=row;
  }
};

template<class MType,class VType>
void vecCmat(VType &v,MType &mat)
{
  //  cout<<v<<endl;
  int row=mat.size(0);
  int col=mat.size(1);
  if(v.size()!=row*col)
  {
    cout<<"Wrong matrix size.Can't covert vector to matrix."<<endl;
    return;
  }
  int pointer=0;

  for(int i=0;i<col;i++)
  {
    mat(LaIndex(),i).inject(v(LaIndex(pointer,pointer+row-1)));
    pointer+=row;
  }

  //cout<<mat.size(0)<<endl;
  //mat.inject(v);
  //cout<<mat.size(0)<<endl;
};


///Add u to the end of v,so that the new v is of size v.size()+u.size()
template<class VType>
void join(VType &v,VType &u)
{
  VType newv(v.size()+u.size());
  //cout<<v.size()+u.size()<<endl;
  if(v.size()!=0)
    newv(LaIndex(0,v.size()-1)).inject(v);
  newv(LaIndex(v.size(),newv.size()-1)).inject(u);
  v=newv;
};


LaGenMatDouble expm(LaGenMatDouble &m);

LaGenMatComplex expm2(LaGenMatDouble &m);

double average(LaGenMatDouble &mat,LaVectorDouble &v);
double average(LaGenMatComplex &mat,LaVectorComplex &v);

COMPLEX average(LaVectorComplex &v1,LaGenMatComplex &mat,LaVectorComplex &v2);

void normalize(LaVectorDouble &v);
void normalize(LaVectorComplex &v);

void normalize(LaGenMatDouble &mat);

COMPLEX Dot_Prod(LaVectorComplex &v1, LaVectorComplex &v2);

//mat+=addmat
void addmat(LaGenMatComplex &mat,const LaGenMatComplex &addmat);
void addmat(LaGenMatDouble &mat,const LaGenMatDouble &addmat);

template<class Type>
void printvector(vector<Type> &v)
{
  int dim=v.size();
  if(dim==0)
    cout<<"Vector size is zero!!!"<<endl;
  else
  {
    for(int i=0;i<dim;i++)
      cout<<v[i]<<" ";
    //cout<<endl;
  }

};

///Read matlab generated operators
template <class MType>
void ReadOperators(ifstream &fin,vector<MType> &OP, int num)
{
	//Initialize OP
    OP.resize(num);
	//cout<<"***********"<<sizeof(typename MType::value_type)<<endl;
	//Read from the file
    for(int i=0;i<num;i++)
    {
		ReadOneOperator(fin,OP[i]);

    }
}

template <class MType>
void ReadOneOperator(ifstream &fin,MType &OP)
{
	int row,col;
	fin.read((char *)&row,sizeof(int));
	fin.read((char *)&col,sizeof(int));
		OP.resize(row,col);
	//cout<<"***********"<<sizeof(typename MType::value_type)<<endl;
	//Read from the file
		fin.read((char *)OP.addr(),sizeof(typename MType::value_type)*row*col);

}

//The following two functions should be added to lacomplex.h
inline bool operator!=(const COMPLEX& _a, const double& _b)
{
  return _a.r != _b || _a.i != 0;
}

inline bool operator==(const COMPLEX& _a, const double& _b)
{
  return _a.r == _b || _a.i == 0;
}

inline COMPLEX  operator+(const COMPLEX& _a, const COMPLEX& _b)
{
  COMPLEX _c;
  _c.r=_a.r+_b.r;
  _c.i=_a.i+_b.i;
  return _c;
}

inline COMPLEX  operator-(const COMPLEX& _a, const COMPLEX& _b)
{
  COMPLEX _c;
  _c.r=_a.r-_b.r;
  _c.i=_a.i-_b.i;
  return _c;
}

inline COMPLEX  operator*(const COMPLEX& _a, const COMPLEX& _b)
{
  COMPLEX _c;
  _c.r=_a.r*_b.r-_a.i*_b.i;
  _c.i=_a.r*_b.i+_a.i*_b.r;
  return _c;
}

inline COMPLEX  operator/(const COMPLEX& _a, const COMPLEX& _b)
{
  COMPLEX _c;
  _c.r=(_a.r*_b.r+_a.i*_b.i)/(_b.i*_b.i+_b.r*_b.r);
  _c.i=(_a.i*_b.r-_a.r*_b.i)/(_b.i*_b.i+_b.r*_b.r);
  return _c;
}


extern "C"
{
  void dsyev_(const char&jobz,const char&uplo,const int&n,double *a,const int&lda,double*w,double *work,int&lwork            ,int&info);
  void zheev_(const char&jobz,const char&uplo,const int&n,COMPLEX *a,const int&lda,double*w,COMPLEX *work,int&lwork,double*rwok,int&info);
  void dgemm_(char *transa, char *transb,long int *m, long int *n, long int *k,
                    double *alpha, const double *a,long int *lda, const double *b,
                   long int *ldb, double *beta, double *c, long int *ldc);
}

void SSMED(double* Matrix,int Dim,double* EigenValue);
void SSMED(COMPLEX* Matrix,int Dim,double* EigenValue);

void blas_mat_mat_mult(double *a,long int arow,long int acol,double *b,long int brow,long int bcol,double *c,long int crow,long int ccol,double alpha=1.0,double beta=0.0);



#endif
