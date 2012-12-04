#include "SuperChain.h"
#include "public.h"


#include <blas1pp.h>
#include <blaspp.h>


snake::physics::SuperChain::SuperChain()
{
  value_type='r';
}


snake::physics::SuperChain::~SuperChain()
{
}


/*!
\fn snake::physics::SuperChain::SuperChain(Block *left,Block *right,Block *oleft,Block *oright,LaGenMatDouble &Hi)
 */
snake::physics::SuperChain::SuperChain(Chain *left,Chain *right,Chain *oleft,Chain *oright,LaGenMatDouble &Hi)
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
    \fn snake::physics::SuperChain::renormwf(DTMat &dtmat)
 */

void snake::physics::SuperChain::renormwfmat(DTMat &dtmat)
{
  //std::cout<<dtmat.trunmat<<std::endl;
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
      std::cout<<"NO HANDSIDE INFORMATION!"<<std::endl;
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
      std::cout<<"NO HANDSIDE INFORMATION!"<<std::endl;
  }
}


/*!
    \fn snake::physics::SuperChain::unrenormwfDTMat &dtmat)
 */
void snake::physics::SuperChain::unrenormwfmat(DTMat &dtmat)
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
    //std::cout<<dtmat.trunmat;
      wfmat=dtmat.trunmat*wfmat;
      rightbase=dtmat.rightbase;
    }
    else
    {
      std::cout<<dtmat.handside<<std::endl;
      std::cout<<"NO HANDSIDE INFORMATION!"<<std::endl;
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
    //std::cout<<dtmat.trunmat;
      wfmatC=dtmat.trunmatC*wfmatC;
      #if TARGET_TWO_WF
      wfmatC2=dtmat.trunmatC*wfmatC2;
      #endif
      rightbase=dtmat.rightbase;
    }
    else
    {
      std::cout<<dtmat.handside<<std::endl;
      std::cout<<"NO HANDSIDE INFORMATION!"<<std::endl;
    }
  }
}





/*!
    \fn snake::physics::SuperChain::prepare()
 */
///This function haven been complished and tested

void snake::physics::SuperChain::prepare()
{
  std::cout<<leftbase<<std::endl;
  std::cout<<rightbase<<std::endl;
  std::cout<<oldleftbase<<std::endl;
  std::cout<<oldrightbase<<std::endl;
  leftbase=L->base;
  rightbase=R->base;
  oldleftbase=oldL->base;
  oldrightbase=oldR->base;
  std::cout<<leftbase<<std::endl;
  std::cout<<rightbase<<std::endl;
  std::cout<<oldleftbase<<std::endl;
  std::cout<<oldrightbase<<std::endl;
  //genindex();
  //genmiddlemap();
}



/*!
    \fn snake::physics::SuperChain::applyop(Rmatrix &op,int thesite)
 */
void snake::physics::SuperChain::applyop(LaGenMatComplex &op,int thesite)
{
  genindex();
    ///Move right
    for(int i=1;i<sitenum-1;i++)
    {
      //std::cout<<"Left site is "<<i<<std::endl;
      //std::cout<<"The mod of wfC is "<<Blas_Norm2(wfC)<<std::endl;
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
      //std::cout<<"The mod of wfC is "<<Blas_Norm2(wfC)<<std::endl;
      moveleft(leftdtmat[i-1],rightdtmat[sitenum-i]);
      oldleftbase=leftdtmat[i-2].tmatbase;
    }
  delindex();
  //fout<<energy;
}


/*!
    \fn snake::physics::SuperChain::write(char *filename)
 */
void snake::physics::SuperChain::write(char *filename)
{
  std::ofstream fout(filename,std::ios_base::out|std::ios_base::binary);
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
    snake::math::writevec(fout,wf);
    wfmat.write(fout);
  }
  else
  {
    snake::math::writevec(fout,wfC);
    wfmatC.write(fout);
    //snake::math::writevec(fout,wfC2);
    //wfmatC2.write(fout);
  }
}


/*!
    \fn snake::physics::SuperChain::read(char *filename)
 */
void snake::physics::SuperChain::read(char *filename)
{
  std::ifstream fin(filename,std::ios_base::in|std::ios_base::binary);
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
    snake::math::readvec(fin,wf);
    wfmat.read(fin);
  }
  else
  {
    snake::math::readvec(fin,wfC);
    wfmatC.read(fin);
    //snake::math::readvec(fin,wfC2);
    //wfmatC2.read(fin);
  }
}


/*!
    \fn snake::physics::SuperChain::applyOPonDot(Rmatrix &OP)
 */
/*
void snake::physics::SuperChain::applyOPonDot(Rmatrix &OP)
{
  normalize(wf);
  evalwfmat(wf,wfmat,TargetGQN);
 // std::cout<<wfmat<<std::endl;
  Rmatrix tempmat(rightbase,leftbase);
  Mat_Mat_Trans_Mult(wfmat,OP,tempmat);
  extractwf(tempmat,wf,tnum2);
  evalwfmat(wf,wfmat,tnum2);

  Rmatrix tempdmat(leftbase,leftbase);
  Rmatrix tempmat2(rightbase,leftbase);
    //std::cout<<wfmatC<<std::endl;
    //std::cout<<freesite.cC[0]<<std::endl;
  Mat_Mat_Mult(wfmat,freesite.n[0],tempmat2);
  Mat_Trans_Mat_Mult(wfmat,tempmat2,tempdmat);
  std::cout<<tempdmat.trace()<<std::endl;
 // std::cout<<wfmat<<std::endl;
}
*/

/*!
\fn snake::physics::SuperChain::renorm()
 */
void snake::physics::SuperChain::renorm(int tn)
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
  //std::cout<<L->dtmat->denmat<<std::endl;

  L->dtmat->findtmat(tn);
  R->dtmat->findtmat(tn);
  L->renorm();
  R->renorm();
  */
}


/*!
\fn snake::physics::SuperChain::renormright()
 */
void snake::physics::SuperChain::renormright(int tn)
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
  std::cout<<"NewRightDim="<<R->base.Dim<<"\t";
}


/*!
\fn snake::physics::SuperChain::renormleft()
 */
void snake::physics::SuperChain::renormleft(int tn)
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
  // std::cout<<L->dtmat->tmatbase<<std::endl;
  // std::cout<<L->dtmat->trunmat<<std::endl;
  L->renorm();
  std::cout<<"NewLeftDim="<<L->base.Dim<<"\t";
}



/*!
\fn snake::physics::SuperChain::calCF()
 */
void snake::physics::SuperChain::calCF(char *filename)
{

  std::ofstream fout(filename);

  evalwfmat(wf,wfmat,TargetGQN);
  renormwfmat(*(L->dtmat));
  renormwfmat(*(R->dtmat));
  wfmat.normalize();
  //std::cout<<wfmat<<std::endl;

  for(int i=L->sitenum-1,j=R->sitenum-1;i>=0&&j>=0;i--,j--)
  {
    L->site[i].eval();
    fout<<corrfunc(L->site[i].c[0],'l',R->site[R->sitenum-1].a[0],'r')<<std::endl;
  }

  fout.close();
}

/*!
\fn snake::physics::SuperChain::corrfunc(LaGenMatDouble &m1,char hand1,LaGenMatDouble &m0,char hand0)
 */
double snake::physics::SuperChain::corrfunc(Rmatrix &m1,char hand1,Rmatrix &m0,char hand0)
{
  //std::cout<<m1<<std::endl;
  //std::cout<<m0<<std::endl;
  Rmatrix tempmat(rightbase,leftbase),tempmat2(rightbase,leftbase);
  ///m0*wf
  //std::cout<<wfmat<<std::endl;
  //std::cout<<m0<<std::endl;
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
 \fn snake::physics::SuperChain::moveright(DTMat &leftdtmat,DTMat &rightdtmat,std::vector<int> &rightordermap)
 */
void snake::physics::SuperChain::moveright(DTMat &leftdtmat,DTMat &rightdtmat)
{
	std::vector<int> temp;

	if(value_type=='r')
	{
		evalwfmat(wf,wfmat,TargetGQN);
		leftdtmat.gendenmat(wfmat,leftbase,rightbase);
		leftdtmat.findtmat(KeptStatesNum);
		leftdtmat.denmat.delmat();
		renormwfmat(leftdtmat);

		LaGenMatDouble wfmatfull;
		wfmat.Convert2Full(wfmatfull);
        snake::math::unordermat(rightbase.ordermap,temp,wfmatfull);
		//std::cout<<wfmat<<std::endl;
		reshapewfmat(wfmatfull,'r');
		//std::cout<<wfmat<<std::endl;
		oldleftbase=leftbase;
		leftbase=kron(oldleftbase,freesite.base);
        snake::math::ordermat(temp,leftbase.ordermap,wfmatfull);
		//std::cout<<wfmatfull<<std::endl;
		wfmat=Rmatrix(wfmatfull,oldrightbase,leftbase,TargetGQN);
		unrenormwfmat(rightdtmat);
		extractwf(wfmat,wf,TargetGQN);
	}
	else
	{
		//std::cout<<wfC.size()<<std::endl;
		evalwfmat(wfC,wfmatC,TargetGQN);
		//std::cout<<wfmatC<<std::endl;
		leftdtmat.gendenmat(wfmatC,leftbase,rightbase);

#if TARGET_TWO_WF
		Cmatrix denmat;
		denmat=leftdtmat.denmatC;
		// std::cout<<denmat<<std::endl;

		evalwfmat(wfC2,wfmatC2,TargetGQN2);
		leftdtmat.gendenmat(wfmatC2,leftbase,rightbase);
		//std::cout<<leftdtmat.denmatC<<std::endl;
		//std::cout<<denmat<<std::endl;
		leftdtmat.denmatC+=denmat;
		//std::cout<<leftdtmat.denmatC<<std::endl;
		//denmat.delmat();
#endif

		leftdtmat.findtmat(KeptStatesNum);
		//std::cout<<leftdtmat.tmatbase.Dim<<" "<<flush;
		//std::cout<<leftdtmat.trunmatC<<std::endl;
		leftdtmat.denmatC.delmat();
		renormwfmat(leftdtmat);
		//std::cout<<wfmatC<<std::endl;

		LaGenMatComplex wfmatfull,gwfmatfull;
		wfmatC.Convert2Full(wfmatfull);
        snake::math::unordermat(rightbase.ordermap,temp,wfmatfull);
		reshapewfmat(wfmatfull,'r');
#if TARGET_TWO_WF
		wfmatC2.Convert2Full(gwfmatfull);
        snake::math::unordermat(rightbase.ordermap,temp,gwfmatfull);
		reshapewfmat(gwfmatfull,'r');
#endif
		//std::cout<<wfmat<<std::endl;
		oldleftbase=leftbase;
		leftbase=kron(oldleftbase,freesite.base);
        snake::math::ordermat(temp,leftbase.ordermap,wfmatfull);
		//std::cout<<wfmatfull<<std::endl;
		//std::cout<<oldrightbase<<std::endl;
		//std::cout<<leftbase<<std::endl;
		wfmatC=Cmatrix(wfmatfull,oldrightbase,leftbase,TargetGQN);
		//std::cout<<wfmatC<<std::endl;
		//wfmatfull.resize(0,0);
#if TARGET_TWO_WF
        snake::math::ordermat(temp,leftbase.ordermap,gwfmatfull);
		wfmatC2=Cmatrix(gwfmatfull,oldrightbase,leftbase,TargetGQN2);
		//gwfmatfull.resize(0,0);
#endif
		unrenormwfmat(rightdtmat);
		//std::cout<<wfmatC<<std::endl;
		extractwf(wfmatC,wfC,TargetGQN);
#if TARGET_TWO_WF
		extractwf(wfmatC2,wfC2,TargetGQN2);
#endif
	}
	//std::cout<<wfmat<<std::endl;
}


/*!
 \fn snake::physics::SuperChain::moveleft(DTMat &leftdtmat,DTMat &rightdtmat,std::vector<int> &rightordermap)
 */
void snake::physics::SuperChain::moveleft(DTMat &leftdtmat,DTMat &rightdtmat)
{
	std::vector<int> temp;

	if(value_type=='r')
	{
		//std::cout<<wf<<std::endl;
		evalwfmat(wf,wfmat,TargetGQN);
		//std::cout<<wfmat<<std::endl;
		rightdtmat.gendenmat(wfmat,leftbase,rightbase);
		//std::cout<<rightdtmat.denmat<<std::endl;
		rightdtmat.findtmat(KeptStatesNum);
		rightdtmat.denmat.delmat();
		renormwfmat(rightdtmat);
		// std::cout<<Mat_Dot_Prod(wfmat,wfmat)<<std::endl;
		LaGenMatDouble wfmatfull;
		wfmat.Convert2Full(wfmatfull);
        snake::math::unordermat(temp,leftbase.ordermap,wfmatfull);
		//printvector(leftbase.ordermap);
		//std::cout<<wfmat<<std::endl;
		//std::cout<<wfmatfull<<std::endl;
		//std::cout<<leftbase<<std::endl;
		//std::cout<<oldleftbase<<std::endl;
		reshapewfmat(wfmatfull,'l');
		//std::cout<<Blas_NormF(wfmatfull)<<std::endl;
		//std::cout<<wfmatfull<<std::endl;
		oldrightbase=rightbase;
		rightbase=kron(freesite.base,oldrightbase);
		//std::cout<<rightbase<<std::endl;
		//std::cout<<oldleftbase<<std::endl;
        snake::math::ordermat(rightbase.ordermap,temp,wfmatfull);
		//std::cout<<wfmatfull<<std::endl;
		//std::cout<<Blas_NormF(wfmatfull)<<std::endl;
		wfmat=Rmatrix(wfmatfull,rightbase,oldleftbase,TargetGQN);
		//std::cout<<wfmat<<std::endl;
		//std::cout<<Mat_Dot_Prod(wfmat,wfmat)<<std::endl;
		unrenormwfmat(leftdtmat);
		//std::cout<<Mat_Dot_Prod(wfmat,wfmat)<<std::endl<<std::endl;
		//std::cout<<wfmat<<std::endl;
		extractwf(wfmat,wf,TargetGQN);
	}
	else
	{
		evalwfmat(wfC,wfmatC,TargetGQN);
		rightdtmat.gendenmat(wfmatC,leftbase,rightbase);
        //std::cout<<wfmatC<<std::endl;
#if TARGET_TWO_WF
		Cmatrix denmat;
		denmat=rightdtmat.denmatC;
		//std::cout<<denmat<<std::endl;
		evalwfmat(wfC2,wfmatC2,TargetGQN2);
		rightdtmat.gendenmat(wfmatC2,leftbase,rightbase);
		rightdtmat.denmatC+=denmat;
		// denmat.delmat();
#endif
		rightdtmat.findtmat(KeptStatesNum);
		//std::cout<<rightdtmat.tmatbase.Dim<<" "<<flush;
		rightdtmat.denmatC.delmat();
		renormwfmat(rightdtmat);
		//std::cout<<wfmatC<<std::endl;
		LaGenMatComplex wfmatfull,gwfmatfull;
		//std::cout<<wfmatC<<std::endl;
		wfmatC.Convert2Full(wfmatfull);
		//std::cout<<wfmatfull<<std::endl;
		//std::cout<<wfmatfull<<std::endl;
        snake::math::unordermat(temp,leftbase.ordermap,wfmatfull);
		reshapewfmat(wfmatfull,'l');
		//std::cout<<wfmatfull<<std::endl;
#if TARGET_TWO_WF
		wfmatC2.Convert2Full(gwfmatfull);
        snake::math::unordermat(temp,leftbase.ordermap,gwfmatfull);
		reshapewfmat(gwfmatfull,'l');
#endif
		//std::cout<<wfmat<<std::endl;
		oldrightbase=rightbase;
		rightbase=kron(freesite.base,oldrightbase);
        snake::math::ordermat(rightbase.ordermap,temp,wfmatfull);
		wfmatC=Cmatrix(wfmatfull,rightbase,oldleftbase,TargetGQN);
		//wfmatfull.resize(0,0);
#if TARGET_TWO_WF
        snake::math::ordermat(rightbase.ordermap,temp,gwfmatfull);
		wfmatC2=Cmatrix(gwfmatfull,rightbase,oldleftbase,TargetGQN2);
		//gwfmatfull.resize(0,0);
#endif
		unrenormwfmat(leftdtmat);
		//std::cout<<wfmatC<<std::endl;
		extractwf(wfmatC,wfC,TargetGQN);
#if TARGET_TWO_WF
		extractwf(wfmatC2,wfC2, TargetGQN2);
#endif
	}

}
