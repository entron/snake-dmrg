#include "public.h"
#include "gqnmat.h"

#include <gmd.h>
#include <laindex.h>
#include <blas1pp.h>
#include <blas2pp.h>
#include <blas3pp.h>
#include <laslv.h>



/**The kronecker product of two LaGenMatDouble
*/
LaGenMatDouble kron(LaGenMatDouble &a,LaGenMatDouble &b)
{
  int arow=a.size(0),brow=b.size(0);
  int acol=a.size(1),bcol=b.size(1);
  LaGenMatDouble c(arow*brow,acol*bcol);
  for(int j=0;j<acol;j++)
    for(int i=0;i<arow;i++)
    {
      int p=i*brow;
      int q=j*bcol;
      for(int n=0;n<bcol;n++)
        for(int m=0;m<brow;m++)
          c(m+p,n+q)=a(i,j)*b(m,n);
    }
  return c;
}


LaGenMatComplex kron(LaGenMatComplex &a,LaGenMatComplex &b)
{
  int arow=a.size(0),brow=b.size(0);
  int acol=a.size(1),bcol=b.size(1);
  LaGenMatComplex c(arow*brow,acol*bcol);
  for(int j=0;j<acol;j++)
    for(int i=0;i<arow;i++)
    {
      int p=i*brow;
      int q=j*bcol;
      for(int n=0;n<bcol;n++)
        for(int m=0;m<brow;m++)
          c(m+p,n+q)=LaComplex(a(i,j))*LaComplex(b(m,n));
    }
  return c;
}

/**The direct sum of two LaGenMatDouble
*/
LaGenMatDouble directsum(LaGenMatDouble &a,LaGenMatDouble &b)
{
  int arow=a.size(0),brow=b.size(0);
  int acol=a.size(1),bcol=b.size(1);
  LaGenMatDouble newmatrix(arow+brow,acol+bcol);
  newmatrix(LaIndex(0,arow+brow-1),LaIndex(0,acol+bcol-1))=0;
  
  if(arow>0&&acol>0)
    newmatrix(LaIndex(0,arow-1),LaIndex(0,acol-1)).inject(a);
  if(brow>0&&bcol>0)
    newmatrix(LaIndex(arow,arow+brow-1),LaIndex(acol,acol+bcol-1)).inject(b); 
  return newmatrix;
}


/**Similar function with chop[] in mathematica
*/
void chop(LaGenMatDouble &a,double err)
{
  int row=a.size(0);
  int col=a.size(1);
  for(int j=0;j<col;j++)
    for(int i=0;i<row;i++)
      if(a(i,j)>-err&&a(i,j)<err)
        a(i,j)=0;
}

void chop(LaGenMatComplex &a,double err)
{
  int row=a.size(0);
  int col=a.size(1);
  double mod;
  for(int j=0;j<col;j++)
    for(int i=0;i<row;i++)
    {
      mod=sqrt(a(i,j).r*a(i,j).r+a(i,j).i*a(i,j).i);
      if(mod<err)
        a(i,j)=LaComplex(0,0);
    }
}




double trace(LaGenMatDouble &mat)
{
  int row=mat.size(0);
  int col=mat.size(1);
  double t=0;
  if(row!=col)
    std::cout<<"Trace only works for square matrices."<<std::endl;
  else
    for(int i=0;i<row;i++)
      t+=mat(i,i);
  return t;
}



/**
*Join an interger and a string togerther as a new string.
*I copy this function from Xiang's program
*/
/*
void combine(char *s1, int ks,char *str) 
{ 
  int km = ks ; 
  int len = strlen( s1 ) + 1  ; 
  while( ks /= 10 ) len++ ; 
  
  str = new char[len + 1] ; 
  str[len] = '\0' ; 
  for ( int i = 0 ; i < strlen(s1) ; i++ )  
    str[i] = s1[i] ; 
  
  do { str[--len] = (km % 10) + '0' ; 
  } while (km /= 10) ; 
  
} 
*/


///Calculate square matrix exponent
LaGenMatDouble expm(LaGenMatDouble &m)
{
  //std::cout<<m<<std::endl;
  int dim=m.size(0);
  LaGenMatDouble r(dim,dim);
  LaGenMatDouble eigvec(dim,dim);
  LaVectorDouble eigval(dim);
  eigvec=m;
  SSMED(eigvec.addr(),dim,eigval.addr());

  for(int i=0;i<dim;i++)
    eigval(i)=exp(eigval(i));
    LaGenMatDouble temp(dim,dim);
  temp=temp.from_diag(eigval);
  //std::cout<<temp<<std::endl;

  Blas_Mat_Mat_Trans_Mult(temp,eigvec,m);
  Blas_Mat_Mat_Mult(eigvec,m,r);
  chop(r,1e-15);
  return r;
}

///Calculate Exp[I*m]
LaGenMatComplex expm2(LaGenMatDouble &m)
{
  //std::cout<<m<<std::endl;
  int dim=m.size(0);
  LaGenMatComplex r(dim,dim);
  LaGenMatDouble eigvec(dim,dim);
  LaGenMatComplex eigvecC(dim,dim);
  LaVectorDouble eigval(dim);
  LaVectorComplex eigvalim(dim);
  eigvec=m;
  SSMED(eigvec.addr(),dim,eigval.addr());

  for(int i=0;i<dim;i++)
    eigvalim(i)=LaComplex(cos(eigval(i)),sin(eigval(i)));
  LaGenMatComplex temp(dim,dim);
  temp=LaComplex(0,0);
  for(int i=0;i<dim;i++)
    temp(i,i)=eigvalim(i);

  chop(temp,1e-15);
  //std::cout<<temp<<std::endl;
  eigvecC=eigvec.to_LaGenMatComplex();
  LaGenMatComplex tempx(dim,dim);
  Blas_Mat_Mat_Trans_Mult(temp,eigvecC,tempx);
  Blas_Mat_Mat_Mult(eigvecC,tempx,r);
  chop(r,1e-15);
  return r;
}

///Solve symmatric matrix eigen problem
void SSMED(double* Matrix,int Dim,double* EigenValue)
{
  assert(Dim>0);
  
  char jobz='V';
  char uplo='U';
  const int n=Dim;
  const int lda=n;
  int info=0;
  
  int lwork=3*Dim;
  
  double*work=new double[lwork];
  assert(work);
  
  dsyev_(jobz,uplo,n,Matrix,lda,EigenValue,work,lwork,info);
  
  delete []work;
  
  //if(info==0) std::cout<<"successful in SSMDiag"<<std::endl;
//
  // else std::cout<<"fail in SSMDiag"<<std::endl;
//
}

void SSMED(COMPLEX* Matrix,int Dim,double* EigenValue)
{
assert(Dim>0);

char jobz='V';
char uplo='U';
const int n=Dim;
const int lda=n;
int info=0;


int lwork=2*Dim;
double*rwork=new double[3*Dim];
assert(rwork);

COMPLEX *work=new COMPLEX[lwork];
assert(work);


zheev_(jobz,uplo,n,Matrix,lda,EigenValue,work,lwork,rwork,info);
delete []rwork;


delete []work;

  //if(info==0) std::cout<<"successful in SSMDiag"<<std::endl;
//
  //else std::cout<<"fail in SSMDiag"<<std::endl;
//
}

double average(LaGenMatDouble &mat,LaVectorDouble &v)
{
  double z;
  z=Blas_Dot_Prod(v,v);
  LaVectorDouble tempv(v);
  Blas_Mat_Vec_Mult(mat,v,tempv);
  return Blas_Dot_Prod(v,tempv)/z;
}

double average(LaGenMatComplex &mat,LaVectorComplex &v)
{
  double z;
  z=Blas_H_Dot_Prod(v,v).r;
  LaVectorComplex tempv(v);
  Blas_Mat_Vec_Mult(mat,v,tempv);
  return Blas_H_Dot_Prod(v,tempv).r/z;
}

///Be carefull that v1,v2 will not be normalized
COMPLEX average(LaVectorComplex &v1,LaGenMatComplex &mat,LaVectorComplex &v2)
{
  COMPLEX z;
  // std::cout<<v2<<std::endl;
  LaVectorComplex tempv(v2.size());
  //std::cout<<tempv<<std::endl;
  Blas_Mat_Vec_Mult(mat,v2,tempv);

  z=Dot_Prod(v1,tempv);
  //std::cout<<LaComplex(z)<<std::endl;
  return z;
}


void normalize(LaVectorDouble &v)
{
  double norm;
  norm=Blas_Norm2(v);
  Blas_Scale(1/norm,v);
}

void normalize(LaVectorComplex &v)
{
  double norm;
  norm=Blas_Norm2(v);
  v.scale(LaComplex(1/norm,0));
}


void normalize(LaGenMatDouble &mat)
{
  double norm;
  norm=Blas_NormF(mat);
  Blas_Scale(1/norm,mat);
}


COMPLEX Dot_Prod(LaVectorComplex &v1, LaVectorComplex &v2)
{
  assert(v1.size()==v2.size());
  COMPLEX result;
  double real=0,imag=0;
  int n=v1.size();
  for(int i=0;i<n;i++)
  {
    real+=(v1(i).r*v2(i).r+v1(i).i*v2(i).i);
    imag+=(v1(i).r*v2(i).i-v1(i).i*v2(i).r);
  }
  result=LaComplex(real,imag);
  return result;
}


void blas_mat_mat_mult(double *a,long int arow,long int acol,double *b,long int brow,long int bcol,double *c,long int crow,long int ccol,double alpha,double beta)
{
  char t='N';
  assert(acol==brow);
  assert(arow==crow);
  assert(bcol==ccol);
  dgemm_(&t, &t, &arow, &bcol, &acol, &alpha, a, &arow, b, 
                 &brow, &beta, c, &crow);
}



void addmat(LaGenMatComplex &mat,const LaGenMatComplex &addmat)
{
  int row=mat.size(0);
  int col=mat.size(1);
  assert(row==addmat.size(0));
  assert(col==addmat.size(1));
  int n=row*col;//
  LaComplex *p,*padd;
  p=(LaComplex*)mat.addr();
  padd=(LaComplex*)addmat.addr();
  for(int i=0;i<n;i++)
  {
    (*p)+=(*padd);
    p++;
    padd++;
  }
}

void addmat(LaGenMatDouble &mat,const LaGenMatDouble &addmat)
{
  int row=mat.size(0);
  int col=mat.size(1);
  assert(row==addmat.size(0));
  assert(col==addmat.size(1));
  int n=row*col;//
  double *p,*padd;
  p=mat.addr();
  padd=addmat.addr();
  for(int i=0;i<n;i++)
  {
    (*p)+=(*padd);
    p++;
    padd++;
  }
}
