#include "supblock.h"

/*!
\fn SupBlock::CalGroundState(GQN tTargetGQN)
 */
void SupBlock::CalGroundState()
{
  //cout<<leftbase<<endl;
  //cout<<rightbase<<endl;
  //cout<<Hlr<<endl;
  //cout<<L->hamiltonian->H<<endl;
  //cout<<R->hamiltonian->H<<endl;
  
  midsite1=allfreesites[L->sitenum-1];
  midsite2=allfreesites[L->sitenum];
  double *eigval,**eigvec;
  int dim=calDim();
  int EigStatesNum=1;

  cout<<"SupL="<<sitenum<<"\t";
  cout<<"TargetGQNNum="; printvector(TargetGQN); cout<<"\t";
  cout<<"SupHamDim="<<dim<<"\t";

  mtimes=0;
  H2Dim=Hlr.size(0);

  eigval = new double[EigStatesNum];
  eigvec = new double*[dim];
  for (int i=0; i<dim; i++) eigvec[i] = new double[EigStatesNum];

  genindex();
  genmiddlemap(TargetGQN);
  dsaupd(dim,EigStatesNum,eigval,eigvec);

  cout<<"MultTimes="<<mtimes<<"\t";
  multnum+=mtimes;
  //cout.width(20);
  cout.precision(15);
  cout<<"Eg="<<eigval[0]<<"\t";

  wf.resize(dim,1);
  for(int i=0;i<dim;i++) wf(i)=eigvec[i][0];

  delete [] eigval;
  for (int i=0; i<dim; i++)
    delete [] eigvec[i];
  delete [] eigvec;
  deletemiddlemap();
  delindex();
}

/*!
\fn SupBlock::calDim()
 */
int SupBlock::calDim()
{
  int Dim=0;
  for(int i=0;i<leftbase.subnum;i++)
    for(int j=0;j<rightbase.subnum;j++)
    {
    	//cout<<leftbase.subgqn[i]+rightbase.subgqn[j]<<endl;
    	//printvector(TargetGQN);
    	if(leftbase.subgqn[i]+rightbase.subgqn[j]==TargetGQN)
    		Dim+=leftbase.dim[i]*rightbase.dim[j];
    }
  return Dim;
}

void SupBlock::av(int n,double *in,double *out)
{
  rightmult(in,out);
  leftmult(in,out);
  middlemult(Hlr,in,out);
  mtimes++;
}


/*!
\fn SupBlock::rightmult(double *in,double *out)
 */
void SupBlock::rightmult(double *in,double *out)
{
  double *instart,*outstart;
  instart=in;
  outstart=out;
  for(int j=0;j<leftbase.subnum;j++)
    for(int i=0;i<rightbase.subnum;i++)
      if(rightbase.subgqn[i]+leftbase.subgqn[j]==TargetGQN)
      {
        blas_mat_mat_mult(R->hamiltonian->H.submat[R->hamiltonian->H.pmat(i,i)].addr(),rightbase.dim[i],rightbase.dim[i],instart,rightbase.dim[i],leftbase.dim[j],outstart,rightbase.dim[i],leftbase.dim[j]);
        instart+=rightbase.dim[i]*leftbase.dim[j];
        outstart+=rightbase.dim[i]*leftbase.dim[j];//Can use instart instead
      }
}


/*!
\fn SupBlock::leftmult(double *in ,double *out)
 */
void SupBlock::leftmult(double *in ,double *out)
{
  double *instart,*outstart;
  instart=in;
  outstart=out;
  for(int j=0;j<leftbase.subnum;j++)
    for(int i=0;i<rightbase.subnum;i++)
      if(rightbase.subgqn[i]+leftbase.subgqn[j]==TargetGQN)
      {
        blas_mat_mat_mult(instart,rightbase.dim[i],leftbase.dim[j],L->hamiltonian->H.submat[L->hamiltonian->H.pmat(j,j)].addr(),leftbase.dim[j],leftbase.dim[j],outstart,rightbase.dim[i],leftbase.dim[j],1.0,1.0);
        instart+=rightbase.dim[i]*leftbase.dim[j];
        outstart+=rightbase.dim[i]*leftbase.dim[j];
      }
}

/*
void SupBlock::middlemult(LaGenMatDouble &TO,double *in,double *out)
{
  for(int i=0;i<H2Dim;i++)
    for(int k=0;k<H2Dim;k++)
      if(TO(i,k)!=0)
        for(int j=0;j<mapdim[i];j++)
          out[map[i][j]]+=TO(i,k)*in[map[k][j]];
}
*/


/*!
\fn SupBlock::genindex()
*/
void SupBlock::genindex()
{
  GQNBase b1=midsite1.base;
  GQNBase b2=midsite2.base;
  //cout<<b<<endl;
  index=new int*** [b1.subnum];
  for(int j=0;j<b1.subnum;j++)
  {
    index[j]=new int** [b1.dim[j]];
    for(int r=0;r<b1.dim[j];r++)
    {
      index[j][r]=new int* [b2.subnum];
      for(int k=0;k<b2.subnum;k++)
        index[j][r][k]=new int [b2.dim[k]];
    }
  }
  int p=0;
  for(int j=0;j<b1.subnum;j++)
    for(int r=0;r<b1.dim[j];r++)
      for(int k=0;k<b2.subnum;k++)
        for(int s=0;s<b2.dim[k];s++)
        {
          index[j][r][k][s]=p;
          p++;
        }
}

/*!
\fn SupBlock::delindex()
*/
void SupBlock::delindex()
{
  GQNBase b1=midsite1.base;
  GQNBase b2=midsite2.base;
  for(int j=0;j<b1.subnum;j++)
  {
    for(int r=0;r<b1.dim[j];r++)
    {
      for(int k=0;k<b2.subnum;k++)
        delete [] index[j][r][k];
      delete [] index[j][r];
    }
    delete [] index[j];
  }
  delete [] index;
}


/*!
\fn SupBlock::reshape(LaGenMatDouble &mat,Site &freesite,char hand)
*/
void SupBlock::genmiddlemap(vector<GQN> &tgqn)
{
  int middleDim=midsite1.base.Dim*midsite2.base.Dim;
  int count[middleDim];
  mapdim=new int [middleDim];
  map=new int *[middleDim];
  for(int i=0;i<middleDim;i++)
    mapdim[i]=0;
  GQNBase leftsite=midsite1.base;
  GQNBase rightsite=midsite2.base;

  for(int m=0;m<leftbase.subnum;m++)
  {
    for(int i=0;i<oldleftbase.subnum;i++)
      for(int p=0;p<oldleftbase.dim[i];p++)
        for(int j=0;j<leftsite.subnum;j++)
          if(leftbase.subgqn[m]==oldleftbase.subgqn[i]+leftsite.subgqn[j])
            for(int r=0;r<leftsite.dim[j];r++)
              for(int k=0;k<rightsite.subnum;k++)
                for(int s=0;s<rightsite.dim[k];s++)
                  for(int l=0;l<oldrightbase.subnum;l++)
                    for(int q=0;q<oldrightbase.dim[l];q++)
                      if(rightsite.subgqn[k]+oldrightbase.subgqn[l]+leftbase.subgqn[m]==tgqn)
                        mapdim[index[j][r][k][s]]+=1;//Exam this by an example then you  will understand it.
  }

  for(int i=0;i<middleDim;i++)
  {
    map[i]=new int [mapdim[i]];
    count[i]=0;
  }

  for(int m=0,pointer=0;m<leftbase.subnum;m++)
  {
    for(int i=0;i<oldleftbase.subnum;i++)
      for(int p=0;p<oldleftbase.dim[i];p++)
        for(int j=0;j<leftsite.subnum;j++)
          if(leftbase.subgqn[m]==oldleftbase.subgqn[i]+leftsite.subgqn[j])
            for(int r=0;r<leftsite.dim[j];r++)
              for(int k=0;k<rightsite.subnum;k++)
                for(int s=0;s<rightsite.dim[k];s++)
                  for(int l=0;l<oldrightbase.subnum;l++)
                    for(int q=0;q<oldrightbase.dim[l];q++)
                      if(rightsite.subgqn[k]+oldrightbase.subgqn[l]+leftbase.subgqn[m]==tgqn)
                      {
                        map[index[j][r][k][s]][count[index[j][r][k][s]]]=pointer;
                        count[index[j][r][k][s]]++;
                        pointer+=1;
                      }
  }
}

/*!
\fn SupBlock::deletemiddlemap()
*/
void SupBlock::deletemiddlemap()
{
  int middleDim=midsite1.base.Dim*midsite2.base.Dim;
  for(int i=0;i<middleDim;i++)
    delete [] map[i];
  delete [] map;
  delete [] mapdim;
}



/*!
\fn SupHam::dsaupd(int n,int nev,double *Evals,double **Evecs)
 */
void SupBlock::dsaupd(int n,int nev,double *Evals,double **Evecs)
{
  int ido = 0;/* Initialization of the reverse communication
		  parameter. */

  char bmat[2] = "I";/* Specifies that the right hand side matrix
			 should be the identity matrix; this makes
			 the problem a standard eigenvalue problem.
			 Setting bmat = "G" would have us solve the
			 problem Av = lBv (this would involve using
			 some other programs from BLAS, however). */
  char which[3] = "SA";/* Ask for the nev eigenvalues of smallest
			   magnitude.  The possible options are
			   LM: largest magnitude
			   SM: smallest magnitude
			   LA: largest real component
			   SA: smallest real compoent
			   LI: largest imaginary component
			   SI: smallest imaginary component */
  double tol = 1e-9;/* Sets the tolerance; tol<=0 specifies
		       machine precision */
  double *resid;
  resid = new double[n];
  int ncv =24;/* The largest number of basis vectors that will
		      be used in the Implicitly Restarted Arnoldi
		      Process.  Work per major iteration is
		      proportional to N*NCV*NCV. */
  if (ncv>n) ncv = n;
  double *v;
  int ldv = n;
  v = new double[ldv*ncv];
  int *iparam;
  iparam = new int[11];/* An array used to pass information to the routines
			   about their functional modes. */
  iparam[0] = 1;// Specifies the shift strategy (1->exact)
  iparam[2] = 3*n;// Maximum number of iterations
  iparam[6] = 1;/* Sets the mode of dsaupd.
		      1 is exact shifting,
		      2 is user-supplied shifts,
		      3 is shift-invert mode,
		      4 is buckling mode,
		      5 is Cayley mode. */

  int *ipntr;
  ipntr = new int[11];/* Indicates the locations in the work array workd
			  where the input and output vectors in the
			  callback routine are located. */
  double *workd;
  workd = new double[3*n];
  double *workl;
  workl = new double[ncv*(ncv+8)];
  int lworkl = ncv*(ncv+8);/* Length of the workl array */
  int info = 0;/* Passes convergence information out of the iteration
		   routine. */
  int rvec = 1;  /* Specifies that eigenvectors should be calculated */
  int *select;
  select = new int[ncv];
  double *d;
  d = new double[2*ncv];/* This vector will return the eigenvalues from
			    the second routine, dseupd. */
  double sigma;
  int ierr;

   /* Here we enter the main loop where the calculations are
     performed.  The communication parameter ido tells us when
     the desired tolerance is reached, and at that point we exit
     and extract the solutions. */

  do {
    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
            &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, &info);

     /* From those results, the eigenvalues and vectors are
     extracted. */
    if ((ido==1)||(ido==-1)) av(n, workd+ipntr[0]-1, workd+ipntr[1]-1);
  } while ((ido==1)||(ido==-1));

  if (info<0) {
    cout << "Error with dsaupd, info = " << info << "\n";
    cout << "Check documentation in dsaupd\n\n";
  } else {
    dseupd_(&rvec, "All", select, d, v, &ldv, &sigma, bmat,
            &n, which, &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, &ierr);

    if (ierr!=0) {
      cout << "Error with dseupd, info = " << ierr << "\n";
      cout << "Check the documentation of dseupd.\n\n";
    } else if (info==1) {
      cout << "Maximum number of iterations reached.\n\n";
    } else if (info==3) {
      cout << "No shifts could be applied during implicit\n";
      cout << "Arnoldi update, try increasing NCV.\n\n";
    }

     /* Before exiting, we copy the solution information over to
       the arrays of the calling program, then clean up the
       memory used by this routine.  For some reason, when I
       don't find the eigenvectors I need to reverse the order of
       the values. */
    int i, j;
    for (i=0; i<nev; i++) Evals[i] = d[i];
    for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[j][i] = v[i*n+j];

    delete [] resid;
    delete [] v;
    delete [] iparam;
    delete [] ipntr;
    delete [] workd;
    delete [] workl;
    delete [] select;
    delete [] d;
  }
}
