#include "supblock.h"
#include "public.h"


#include <blas1pp.h>
#include <blaspp.h>


SupBlock::SupBlock()
{
  value_type='r';
}


SupBlock::~SupBlock()
{
}


/*!
\fn SupBlock::SupBlock(Block *left,Block *right,Block *oleft,Block *oright,LaGenMatDouble &Hi)
 */
SupBlock::SupBlock(Block *left,Block *right,Block *oleft,Block *oright,LaGenMatDouble &Hi)
{
  value_type='r';
  L=left;
  R=right;
  oldL=oleft;
  oldR=oright;
  leftbase=L->base;
  rightbase=R->base;
  oldleftbase=oldL->base;
  oldrightbase=oldR->base;
  sitenum=L->sitenum+R->sitenum;
  Hlr=Hi;
}


/*!
    \fn SupBlock::renormwf(DTMat &dtmat)
 */

void SupBlock::renormwfmat(DTMat &dtmat)
{
  //cout<<dtmat.trunmat<<endl;
  if(value_type=='r')
  {
    if(dtmat.handside=='l')
    {
      wfmat=wfmat*dtmat.trunmat;
      leftbase=dtmat.tmatbase;
    }
    else if(dtmat.handside=='r')
    {
        Rmatrix wfmattemp=wfmat;
        wfmat.resize(dtmat.tmatbase,wfmat.colbase);
        Mat_Trans_Mat_Mult(dtmat.trunmat,wfmattemp,wfmat);
      //wfmat=dtmat.trunmattrans*wfmat;
      rightbase=dtmat.tmatbase;
    }
    else
      cout<<"NO HANDSIDE INFORMATION!"<<endl;
  }
  else
  {
    if(dtmat.handside=='l')
    {
      wfmatC=wfmatC*dtmat.trunmatC;
      #if TARGET_TWO_WF
      wfmatC2=wfmatC2*dtmat.trunmatC;
      #endif
      leftbase=dtmat.tmatbase;
    }
    else if(dtmat.handside=='r')
    {
        Cmatrix wfmattemp=wfmatC;
        wfmatC.resize(dtmat.tmatbase,wfmatC.colbase);
        Mat_Trans_Mat_Mult(dtmat.trunmatC,wfmattemp,wfmatC);
      //wfmatC=dtmat.trunmattrans2*wfmatC;
      #if TARGET_TWO_WF
      wfmattemp=wfmatC2;
      wfmatC2.resize(dtmat.tmatbase,wfmatC2.colbase);
      Mat_Trans_Mat_Mult(dtmat.trunmatC,wfmattemp,wfmatC2);
      #endif
      rightbase=dtmat.tmatbase;
    }
    else
      cout<<"NO HANDSIDE INFORMATION!"<<endl;
  }
}


/*!
    \fn SupBlock::unrenormwfDTMat &dtmat)
 */
void SupBlock::unrenormwfmat(DTMat &dtmat)
{
  if(value_type=='r')
  {
    if(dtmat.handside=='l')
    {
      Rmatrix wfmattemp=wfmat;
      wfmat.resize(wfmat.rowbase, dtmat.trunmat.rowbase);
      Mat_Mat_Trans_Mult(wfmattemp,dtmat.trunmat,wfmat);
      //wfmat=wfmat*dtmat.trunmattrans;
      leftbase=dtmat.leftbase;
    }
    else if(dtmat.handside=='r')
    {
    //cout<<dtmat.trunmat;
      wfmat=dtmat.trunmat*wfmat;
      rightbase=dtmat.rightbase;
    }
    else
    {
      cout<<dtmat.handside<<endl;
      cout<<"NO HANDSIDE INFORMATION!"<<endl;
    }
  }
  else
  {
    if(dtmat.handside=='l')
    {
        Cmatrix wfmattemp=wfmatC;
        wfmatC.resize(wfmatC.rowbase, dtmat.trunmatC.rowbase);
        Mat_Mat_Trans_Mult(wfmattemp,dtmat.trunmatC,wfmatC);
      //wfmatC=wfmatC*dtmat.trunmattrans2;
      #if TARGET_TWO_WF
      wfmattemp=wfmatC2;
      wfmatC2.resize(wfmatC2.rowbase, dtmat.trunmatC.rowbase);
      Mat_Mat_Trans_Mult(wfmattemp,dtmat.trunmatC,wfmatC2);
      #endif
      leftbase=dtmat.leftbase;
    }
    else if(dtmat.handside=='r')
    {
    //cout<<dtmat.trunmat;
      wfmatC=dtmat.trunmatC*wfmatC;
      #if TARGET_TWO_WF
      wfmatC2=dtmat.trunmatC*wfmatC2;
      #endif
      rightbase=dtmat.rightbase;
    }
    else
    {
      cout<<dtmat.handside<<endl;
      cout<<"NO HANDSIDE INFORMATION!"<<endl;
    }
  }
}





/*!
    \fn SupBlock::prepare()
 */
///This function haven been complished and tested

void SupBlock::prepare()
{
  cout<<leftbase<<endl;
  cout<<rightbase<<endl;
  cout<<oldleftbase<<endl;
  cout<<oldrightbase<<endl;
  leftbase=L->base;
  rightbase=R->base;
  oldleftbase=oldL->base;
  oldrightbase=oldR->base;
  cout<<leftbase<<endl;
  cout<<rightbase<<endl;
  cout<<oldleftbase<<endl;
  cout<<oldrightbase<<endl;
  //genindex();
  //genmiddlemap();
}



/*!
    \fn SupBlock::applyop(Rmatrix &op,int thesite)
 */
void SupBlock::applyop(LaGenMatComplex &op,int thesite)
{
  genindex();
    ///Move right
    for(int i=1;i<sitenum-1;i++)
    {
      //cout<<"Left site is "<<i<<endl;
      //cout<<"The mod of wfC is "<<Blas_Norm2(wfC)<<endl;
      if(i==thesite)
      {
      genmiddlemap(TargetGQN);
      middlemult(op,wfC);
      }
      moveright(leftdtmat[i],rightdtmat[sitenum-i-1]);
      oldrightbase=rightdtmat[sitenum-i-2].tmatbase;
    }

  if(thesite==sitenum-1)
  {
    genmiddlemap(TargetGQN);
    middlemult(op,wfC);
  }

    ///Move left
    for(int i=sitenum-1;i>1;i--)
    {
      //cout<<"The mod of wfC is "<<Blas_Norm2(wfC)<<endl;
      moveleft(leftdtmat[i-1],rightdtmat[sitenum-i]);
      oldleftbase=leftdtmat[i-2].tmatbase;
    }
  delindex();
  //fout<<energy;
}


/*!
    \fn SupBlock::write(char *filename)
 */
void SupBlock::write(char *filename)
{
  ofstream fout(filename,ios_base::out|ios_base::binary);
  fout.write((char*)&sitenum,sizeof sitenum);
  int TargetGQNNum;
  TargetGQNNum=TargetGQN.size();
  fout.write((char*)&TargetGQNNum,sizeof TargetGQNNum);
  for(int i=0;i<TargetGQNNum;i++)
	  TargetGQN[i].write(fout);
  fout.write((char*)&KeptStatesNum,sizeof KeptStatesNum);
  fout.write(&value_type,sizeof value_type);
  rightbase.write(fout);
  leftbase.write(fout);
  for(int i=0;i<sitenum;i++)
  {
    leftdtmat[i].write(fout);
    rightdtmat[i].write(fout);
  }
  if(value_type=='r')
  {
    writevec(fout,wf);
    wfmat.write(fout);
  }
  else
  {
    writevec(fout,wfC);
    wfmatC.write(fout);
    //writevec(fout,wfC2);
    //wfmatC2.write(fout);
  }
}


/*!
    \fn SupBlock::read(char *filename)
 */
void SupBlock::read(char *filename)
{
  ifstream fin(filename,ios_base::in|ios_base::binary);
  fin.read((char*)&sitenum,sizeof sitenum);
  int TargetGQNNum;
  fin.read((char*)&TargetGQNNum,sizeof(int));
  TargetGQN.resize(TargetGQNNum);
  for(int i=0;i<TargetGQNNum;i++)
	  TargetGQN[i].read(fin);

  fin.read((char*)&KeptStatesNum,sizeof KeptStatesNum);
  fin.read(&value_type,sizeof value_type);
  rightbase.read(fin);
  leftbase.read(fin);
  for(int i=0;i<sitenum;i++)
  {
    leftdtmat[i].read(fin);
    rightdtmat[i].read(fin);
  }
  if(value_type=='r')
  {
    readvec(fin,wf);
    wfmat.read(fin);
  }
  else
  {
    readvec(fin,wfC);
    wfmatC.read(fin);
    //readvec(fin,wfC2);
    //wfmatC2.read(fin);
  }
}


/*!
    \fn SupBlock::applyOPonDot(Rmatrix &OP)
 */
/*
void SupBlock::applyOPonDot(Rmatrix &OP)
{
  normalize(wf);
  evalwfmat(wf,wfmat,TargetGQN);
 // cout<<wfmat<<endl;
  Rmatrix tempmat(rightbase,leftbase);
  Mat_Mat_Trans_Mult(wfmat,OP,tempmat);
  extractwf(tempmat,wf,tnum2);
  evalwfmat(wf,wfmat,tnum2);

  Rmatrix tempdmat(leftbase,leftbase);
  Rmatrix tempmat2(rightbase,leftbase);
    //cout<<wfmatC<<endl;
    //cout<<freesite.cC[0]<<endl;
  Mat_Mat_Mult(wfmat,freesite.n[0],tempmat2);
  Mat_Trans_Mat_Mult(wfmat,tempmat2,tempdmat);
  cout<<tempdmat.trace()<<endl;
 // cout<<wfmat<<endl;
}
*/

/*!
\fn SupBlock::renorm()
 */
void SupBlock::renorm(int tn)
{
renormleft(tn);
renormright(tn);

	/*Use the following code to increase some speed.
  if(value_type=='r')
  {
    evalwfmat(wf,wfmat,TargetGQN);
    L->dtmat->gendenmat(wfmat,L->base,R->base);
    R->dtmat->gendenmat(wfmat,L->base,R->base);
  }
  else
  {
    evalwfmat(wfC,wfmatC,TargetGQN);
    L->dtmat->gendenmat(wfmatC,L->base,R->base);
    R->dtmat->gendenmat(wfmatC,L->base,R->base);
  }
  //cout<<L->dtmat->denmat<<endl;

  L->dtmat->findtmat(tn);
  R->dtmat->findtmat(tn);
  L->renorm();
  R->renorm();
  */
}


/*!
\fn SupBlock::renormright()
 */
void SupBlock::renormright(int tn)
{
  if(value_type=='r')
  {
    evalwfmat(wf,wfmat,TargetGQN);
    R->dtmat->gendenmat(wfmat,L->base,R->base);
  }
  else
  {
    evalwfmat(wfC,wfmatC,TargetGQN);
    R->dtmat->gendenmat(wfmatC,L->base,R->base);
  }
  R->dtmat->findtmat(tn);
  R->renorm();
  cout<<"NewRightDim="<<R->base.Dim<<"\t";
}


/*!
\fn SupBlock::renormleft()
 */
void SupBlock::renormleft(int tn)
{
  if(value_type=='r')
  {
    evalwfmat(wf,wfmat,TargetGQN);
    L->dtmat->gendenmat(wfmat,L->base,R->base);
  }
  else
  {
    evalwfmat(wfC,wfmatC,TargetGQN);
    L->dtmat->gendenmat(wfmatC,L->base,R->base);
  }

  L->dtmat->findtmat(tn);
  // cout<<L->dtmat->tmatbase<<endl;
  // cout<<L->dtmat->trunmat<<endl;
  L->renorm();
  cout<<"NewLeftDim="<<L->base.Dim<<"\t";
}



/*!
\fn SupBlock::calCF()
 */
void SupBlock::calCF(char *filename)
{

  ofstream fout(filename);

  evalwfmat(wf,wfmat,TargetGQN);
  renormwfmat(*(L->dtmat));
  renormwfmat(*(R->dtmat));
  wfmat.normalize();
  //cout<<wfmat<<endl;

  for(int i=L->sitenum-1,j=R->sitenum-1;i>=0&&j>=0;i--,j--)
  {
    L->site[i].eval();
    fout<<corrfunc(L->site[i].c[0],'l',R->site[R->sitenum-1].a[0],'r')<<endl;
  }

  fout.close();
}

/*!
\fn SupBlock::corrfunc(LaGenMatDouble &m1,char hand1,LaGenMatDouble &m0,char hand0)
 */
double SupBlock::corrfunc(Rmatrix &m1,char hand1,Rmatrix &m0,char hand0)
{
  //cout<<m1<<endl;
  //cout<<m0<<endl;
  Rmatrix tempmat(rightbase,leftbase),tempmat2(rightbase,leftbase);
  ///m0*wf
  //cout<<wfmat<<endl;
  //cout<<m0<<endl;
  if(hand0=='l')
    Mat_Mat_Mult(wfmat,m0,tempmat);
  else
    Mat_Trans_Mat_Mult(m0,wfmat,tempmat);

  ///m1*(m0*wf)
  if(hand1=='l')
    Mat_Mat_Mult(tempmat,m1,tempmat2);
  else
    Mat_Trans_Mat_Mult(m1,tempmat,tempmat2);

  ///Calculate corr.Note that wf should be normalized before evalwfmat()
  return Mat_Dot_Prod(wfmat,tempmat2);
}

/*!
 \fn SupBlock::moveright(DTMat &leftdtmat,DTMat &rightdtmat,vector<int> &rightordermap)
 */
void SupBlock::moveright(DTMat &leftdtmat,DTMat &rightdtmat)
{
	vector<int> temp;

	if(value_type=='r')
	{
		evalwfmat(wf,wfmat,TargetGQN);
		leftdtmat.gendenmat(wfmat,leftbase,rightbase);
		leftdtmat.findtmat(KeptStatesNum);
		leftdtmat.denmat.delmat();
		renormwfmat(leftdtmat);

		LaGenMatDouble wfmatfull;
		wfmat.Convert2Full(wfmatfull);
		unordermat(rightbase.ordermap,temp,wfmatfull);
		//cout<<wfmat<<endl;
		reshapewfmat(wfmatfull,'r');
		//cout<<wfmat<<endl;
		oldleftbase=leftbase;
		leftbase=kron(oldleftbase,freesite.base);
		ordermat(temp,leftbase.ordermap,wfmatfull);
		//cout<<wfmatfull<<endl;
		wfmat=Rmatrix(wfmatfull,oldrightbase,leftbase,TargetGQN);
		unrenormwfmat(rightdtmat);
		extractwf(wfmat,wf,TargetGQN);
	}
	else
	{
		//cout<<wfC.size()<<endl;
		evalwfmat(wfC,wfmatC,TargetGQN);
		//cout<<wfmatC<<endl;
		leftdtmat.gendenmat(wfmatC,leftbase,rightbase);

#if TARGET_TWO_WF
		Cmatrix denmat;
		denmat=leftdtmat.denmatC;
		// cout<<denmat<<endl;

		evalwfmat(wfC2,wfmatC2,TargetGQN2);
		leftdtmat.gendenmat(wfmatC2,leftbase,rightbase);
		//cout<<leftdtmat.denmatC<<endl;
		//cout<<denmat<<endl;
		leftdtmat.denmatC+=denmat;
		//cout<<leftdtmat.denmatC<<endl;
		//denmat.delmat();
#endif

		leftdtmat.findtmat(KeptStatesNum);
		//cout<<leftdtmat.tmatbase.Dim<<" "<<flush;
		//cout<<leftdtmat.trunmatC<<endl;
		leftdtmat.denmatC.delmat();
		renormwfmat(leftdtmat);
		//cout<<wfmatC<<endl;

		LaGenMatComplex wfmatfull,gwfmatfull;
		wfmatC.Convert2Full(wfmatfull);
		unordermat(rightbase.ordermap,temp,wfmatfull);
		reshapewfmat(wfmatfull,'r');
#if TARGET_TWO_WF
		wfmatC2.Convert2Full(gwfmatfull);
		unordermat(rightbase.ordermap,temp,gwfmatfull);
		reshapewfmat(gwfmatfull,'r');
#endif
		//cout<<wfmat<<endl;
		oldleftbase=leftbase;
		leftbase=kron(oldleftbase,freesite.base);
		ordermat(temp,leftbase.ordermap,wfmatfull);
		//cout<<wfmatfull<<endl;
		//cout<<oldrightbase<<endl;
		//cout<<leftbase<<endl;
		wfmatC=Cmatrix(wfmatfull,oldrightbase,leftbase,TargetGQN);
		//cout<<wfmatC<<endl;
		//wfmatfull.resize(0,0);
#if TARGET_TWO_WF
		ordermat(temp,leftbase.ordermap,gwfmatfull);
		wfmatC2=Cmatrix(gwfmatfull,oldrightbase,leftbase,TargetGQN2);
		//gwfmatfull.resize(0,0);
#endif
		unrenormwfmat(rightdtmat);
		//cout<<wfmatC<<endl;
		extractwf(wfmatC,wfC,TargetGQN);
#if TARGET_TWO_WF
		extractwf(wfmatC2,wfC2,TargetGQN2);
#endif
	}
	//cout<<wfmat<<endl;
}


/*!
 \fn SupBlock::moveleft(DTMat &leftdtmat,DTMat &rightdtmat,vector<int> &rightordermap)
 */
void SupBlock::moveleft(DTMat &leftdtmat,DTMat &rightdtmat)
{
	vector<int> temp;

	if(value_type=='r')
	{
		//cout<<wf<<endl;
		evalwfmat(wf,wfmat,TargetGQN);
		//cout<<wfmat<<endl;
		rightdtmat.gendenmat(wfmat,leftbase,rightbase);
		//cout<<rightdtmat.denmat<<endl;
		rightdtmat.findtmat(KeptStatesNum);
		rightdtmat.denmat.delmat();
		renormwfmat(rightdtmat);
		// cout<<Mat_Dot_Prod(wfmat,wfmat)<<endl;
		LaGenMatDouble wfmatfull;
		wfmat.Convert2Full(wfmatfull);
		unordermat(temp,leftbase.ordermap,wfmatfull);
		//printvector(leftbase.ordermap);
		//cout<<wfmat<<endl;
		//cout<<wfmatfull<<endl;
		//cout<<leftbase<<endl;
		//cout<<oldleftbase<<endl;
		reshapewfmat(wfmatfull,'l');
		//cout<<Blas_NormF(wfmatfull)<<endl;
		//cout<<wfmatfull<<endl;
		oldrightbase=rightbase;
		rightbase=kron(freesite.base,oldrightbase);
		//cout<<rightbase<<endl;
		//cout<<oldleftbase<<endl;
		ordermat(rightbase.ordermap,temp,wfmatfull);
		//cout<<wfmatfull<<endl;
		//cout<<Blas_NormF(wfmatfull)<<endl;
		wfmat=Rmatrix(wfmatfull,rightbase,oldleftbase,TargetGQN);
		//cout<<wfmat<<endl;
		//cout<<Mat_Dot_Prod(wfmat,wfmat)<<endl;
		unrenormwfmat(leftdtmat);
		//cout<<Mat_Dot_Prod(wfmat,wfmat)<<endl<<endl;
		//cout<<wfmat<<endl;
		extractwf(wfmat,wf,TargetGQN);
	}
	else
	{
		evalwfmat(wfC,wfmatC,TargetGQN);
		rightdtmat.gendenmat(wfmatC,leftbase,rightbase);
        //cout<<wfmatC<<endl;
#if TARGET_TWO_WF
		Cmatrix denmat;
		denmat=rightdtmat.denmatC;
		//cout<<denmat<<endl;
		evalwfmat(wfC2,wfmatC2,TargetGQN2);
		rightdtmat.gendenmat(wfmatC2,leftbase,rightbase);
		rightdtmat.denmatC+=denmat;
		// denmat.delmat();
#endif
		rightdtmat.findtmat(KeptStatesNum);
		//cout<<rightdtmat.tmatbase.Dim<<" "<<flush;
		rightdtmat.denmatC.delmat();
		renormwfmat(rightdtmat);
		//cout<<wfmatC<<endl;
		LaGenMatComplex wfmatfull,gwfmatfull;
		//cout<<wfmatC<<endl;
		wfmatC.Convert2Full(wfmatfull);
		//cout<<wfmatfull<<endl;
		//cout<<wfmatfull<<endl;
		unordermat(temp,leftbase.ordermap,wfmatfull);
		reshapewfmat(wfmatfull,'l');
		//cout<<wfmatfull<<endl;
#if TARGET_TWO_WF
		wfmatC2.Convert2Full(gwfmatfull);
		unordermat(temp,leftbase.ordermap,gwfmatfull);
		reshapewfmat(gwfmatfull,'l');
#endif
		//cout<<wfmat<<endl;
		oldrightbase=rightbase;
		rightbase=kron(freesite.base,oldrightbase);
		ordermat(rightbase.ordermap,temp,wfmatfull);
		wfmatC=Cmatrix(wfmatfull,rightbase,oldleftbase,TargetGQN);
		//wfmatfull.resize(0,0);
#if TARGET_TWO_WF
		ordermat(rightbase.ordermap,temp,gwfmatfull);
		wfmatC2=Cmatrix(gwfmatfull,rightbase,oldleftbase,TargetGQN2);
		//gwfmatfull.resize(0,0);
#endif
		unrenormwfmat(leftdtmat);
		//cout<<wfmatC<<endl;
		extractwf(wfmatC,wfC,TargetGQN);
#if TARGET_TWO_WF
		extractwf(wfmatC2,wfC2, TargetGQN2);
#endif
	}

}
