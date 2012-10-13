#include "supblock.h"


/*!
 \fn SupBlock::evolve(vector<LaGenMatComplex> &rt_td_impurity_OP,vector<LaGenMatComplex> &rt_OP, int timesteps)
 */
void SupBlock::evolve(vector<LaGenMatComplex> &rt_td_impurity_OP,vector<LaGenMatComplex> &rt_OP, int timesteps)
{
	/////////////////////Output reduced density matrix////////////////////////
	string filename;
	filename="./results/rdm.dat";
	ofstream fout_rdm(filename.c_str());
	fout_rdm.width(20);
	fout_rdm.precision(15);
	//////////////////////////////////////////////////////////////////////////
	
	COMPLEX ndot;
	
  //genindex();
	//cout<<wfC<<endl;
	//cout<<wfmatC<<endl;
	normalize(wfC);
	normalize(wfC2);
	//cout<<wfC<<endl;
	
	double last_step_error=totaltrunerror;
	///Time evolve
	for(int n=0;n<timesteps;n++)
	{
		cout<<"TimeStep="<<n<<endl;
		onetimestep(rt_td_impurity_OP[n],rt_OP);
		///vonNeumann entropy when there are two sites on the left
		for(int i=2;i<sitenum-2;i++)
		{
			fout_entropy_t<<rightdtmat[sitenum-i].entropy<<" ";
		}
		fout_entropy_t<<endl;
		
		evalwfmat(wfC,wfmatC,TargetGQN);
#if TARGET_TWO_WF
		evalwfmat(wfC2,wfmatC2,TargetGQN2);
#endif
		//cout<<wfmatC<<endl; cout<<wfmatC2<<endl;
		
		Cmatrix ReducedDM(leftbase,leftbase);
		Mat_Trans_Mat_Mult(wfmatC, wfmatC, ReducedDM);
		//cout<<ReducedDM<<endl;
#if TARGET_TWO_WF
		
		Cmatrix ReducedDM2(leftbase,leftbase);
		Mat_Trans_Mat_Mult(wfmatC2, wfmatC2, ReducedDM2);
		ReducedDM+=ReducedDM2;
		//cout<<ReducedDM2<<endl;
		//Cmatrix ReducedDM3(leftbase,leftbase);
		//cout<<wfmatC2<<endl;
		//cout<<wfmatC<<endl;
		//cout<<ReducedDM3<<endl;
		ReducedDM2=0;
		Mat_Trans_Mat_Mult(wfmatC2, wfmatC, ReducedDM2);
		//cout<<ReducedDM2<<endl;
		ReducedDM+=ReducedDM2;
		
		ReducedDM2=0;
		Mat_Trans_Mat_Mult(wfmatC, wfmatC2, ReducedDM2);
		ReducedDM+=ReducedDM2;
		//cout<<ReducedDM<<endl;
		//fout_rdm<<ReducedDM.submat[ReducedDM.pmat(0,0)](0,0)<<" ";
		//fout_rdm<<ReducedDM.submat[ReducedDM.pmat(1,0)](0,0)<<" ";
		//fout_rdm<<ReducedDM.submat[ReducedDM.pmat(0,1)](0,0)<<" ";
		//fout_rdm<<ReducedDM.submat[ReducedDM.pmat(1,1)](0,0)<<endl;
		
		fout_rdm<<ReducedDM.submat[ReducedDM.pmat(1,0)](0,0).r<<" ";
		fout_rdm<<ReducedDM.submat[ReducedDM.pmat(1,0)](0,0).i<<endl;
		
		
#endif
  	//cout<<ReducedDM<<endl;
		ndot=ReducedDM.submat[ReducedDM.pmat(0,0)](0,0);
	 	COMPLEX offdiag1=ReducedDM.submat[ReducedDM.pmat(0,0)](0,1);
		COMPLEX offdiag2=ReducedDM.submat[ReducedDM.pmat(0,0)](1,0);
		//cout<<"offdiag1="<<offdiag1<<"\t|"<<"offdiag2="<<offdiag2<<"\t";
		double sigmaz=1-2*ndot.r;
		fout_1stsite_n_t<<sigmaz<<endl;
		cout<<"sigma_z(t)="<<sigmaz<<"\t";
		fout_rdm<<offdiag1.r<<" "<<offdiag1.i<<endl;
		cout<<"TotalError="<<totaltrunerror<<"\t";
		fout_steperror_t<<totaltrunerror-last_step_error<<endl;
		last_step_error=totaltrunerror;
		cout<<endl<<endl;
		
	}
  //delindex();
}

void SupBlock::onetimestep(LaGenMatComplex& impurity_OP_t, vector<LaGenMatComplex>& rt_OP)
{
	///Move right
	//cout<<"LKeptStatesNum="<<flush;
	for(int i=1;i<=sitenum-1;i++)
	{
		if(i==1)
		{
  		//cout<<impurity_OP_t<<endl;
  		//cout<<TargetGQN[0]<<endl;
  		midsite1=allfreesites[0];
  		midsite2=allfreesites[1];
  		        genindex();
			genmiddlemap(TargetGQN);
			middlemult(impurity_OP_t,wfC);
			deletemiddlemap();
#if TARGET_TWO_WF

			genmiddlemap(TargetGQN2);
			middlemult(impurity_OP_t,wfC2);
			deletemiddlemap();
#endif
  		        delindex();
		}
		else
		{
  		midsite1=allfreesites[i-1];
  		midsite2=allfreesites[i];
  		genindex();
			genmiddlemap(TargetGQN);
			middlemult(rt_OP[i-1],wfC);
			deletemiddlemap();
#if TARGET_TWO_WF
			genmiddlemap(TargetGQN2);
			middlemult(rt_OP[i-1],wfC2);
			deletemiddlemap();
  		delindex();
#endif
		}
		if(i<=sitenum-2)
		{
			freesite=allfreesites[i+1];
			moveright(leftdtmat[i],rightdtmat[sitenum-i-1]);
			oldrightbase=rightdtmat[sitenum-i-2].tmatbase;
		}
		
	}
	cout<<endl;
	
	
	///Move left
	//cout<<"RKeptStatesNum="<<flush;
	for(int i=sitenum-1;i>=1;i--)
	{
		if(i==1)
		{
  		//cout<<TargetGQN[0]<<endl;
  		midsite1=allfreesites[0];
  		midsite2=allfreesites[1];
  		genindex();
			genmiddlemap(TargetGQN);
			middlemult(impurity_OP_t,wfC);
			deletemiddlemap();
#if TARGET_TWO_WF
			genmiddlemap(TargetGQN2);
			middlemult(impurity_OP_t,wfC2);
			deletemiddlemap();
#endif
  		delindex();
		}
		else
		{
  		midsite1=allfreesites[i-1];
  		midsite2=allfreesites[i];
  		genindex();
			genmiddlemap(TargetGQN);
			middlemult(rt_OP[i-1],wfC);
			deletemiddlemap();
#if TARGET_TWO_WF
			genmiddlemap(TargetGQN2);
			middlemult(rt_OP[i-1],wfC2);
			deletemiddlemap();
#endif
  		delindex();
		}
		
		if(i>=2)
		{
			freesite=allfreesites[i];
			moveleft(leftdtmat[i-1],rightdtmat[sitenum-i]);
			oldleftbase=leftdtmat[i-2].tmatbase;
		}
		
	}
	cout<<endl;
	normalize(wfC);
	normalize(wfC2);
}





void SupBlock::creatoutputfiles()
{
	string filename;
#if CAL_DURING_TIME_EVOLVE
	filename="./results/sigmaz_t.dat";
	fout_1stsite_n_t.open(filename.c_str());
	fout_1stsite_n_t.width(20);
	fout_1stsite_n_t.precision(15);
#endif
	filename="./results/vonneumannentropy_t.dat";
	fout_entropy_t.open(filename.c_str());
	fout_entropy_t.width(20);
	fout_entropy_t.precision(15);
	
	filename="./results/steperror_t.dat";
	fout_steperror_t.open(filename.c_str());
	fout_steperror_t.width(20);
	fout_steperror_t.precision(15);
}

void SupBlock::closeoutputfiles()
{
#if CAL_DURING_TIME_EVOLVE
	fout_1stsite_n_t.close();
#endif
	fout_entropy_t.close();
	fout_steperror_t.close();
	
}

/*!
 \fn SupBlock::loaddtmats(int n)
 */
void SupBlock::loaddtmats()
{
	rightdtmat.resize(sitenum);
	leftdtmat.resize(sitenum);
	
	
	for(int i=1;i<=sitenum/2;i++)
	{
		string filename;
		stringstream stemp2;
		stemp2<<i;
		filename="./data/L";
		filename+=stemp2.str();
		Block lb(filename);
		filename="./data/R";
		filename+=stemp2.str();
		Block rb(filename);
		leftdtmat[i]=*lb.dtmat;
		rightdtmat[i]=*rb.dtmat;
		if(i==1)
		{
			leftdtmat[i].tmatbase=lb.base;
			leftdtmat[i].leftbase=lb.base;
			leftdtmat[i].trunmat.geneye(lb.base);
			//leftdtmat[i].trunmattrans.geneye(lb.base);
			rightdtmat[i].tmatbase=rb.base;
			rightdtmat[i].rightbase=rb.base;
			rightdtmat[i].trunmat.geneye(rb.base);
			//rightdtmat[i].trunmattrans.geneye(rb.base);
		}
	}
	
	
	//cout<<leftdtmat[1].tmatbase<<endl;
	rightdtmat[0].tmatbase.genvacuumbase();
	rightdtmat[0].rightbase.genvacuumbase();
	leftdtmat[0].leftbase.genvacuumbase();
	leftdtmat[0].tmatbase=rightdtmat[0].tmatbase;
	leftdtmat[0].trunmat.geneye(leftdtmat[0].tmatbase);
	//leftdtmat[0].trunmattrans.geneye(leftdtmat[0].tmatbase);
	rightdtmat[0].trunmat.geneye(leftdtmat[0].tmatbase);
	//rightdtmat[0].trunmattrans.geneye(leftdtmat[0].tmatbase);
	
	for(int i=0;i<sitenum;i++)
	{
		rightdtmat[i].handside='r';
		leftdtmat[i].handside='l';
	}
}

/*!
 \fn SupBlock::sweep2left()
 */
void SupBlock::sweep2left(int NewLeftChainLength)
{
	for(int i=NewLeftChainLength;i>1;i--)
	{
		//cout<<"Left site is "<<i<<endl;
		//cout<<wfmat<<endl;
		//cout<<"Norm="<<Blas_Norm2(wf)<<endl;
  	        freesite=allfreesites[i];
		moveleft(leftdtmat[i-1],rightdtmat[sitenum-i]);
		oldleftbase=leftdtmat[i-2].tmatbase;
		
	}
}

void SupBlock::sweep2leftmost()
{
	
	//	cout<<leftbase<<endl;
	//	cout<<rightbase<<endl;
	//	cout<<oldleftbase<<endl;
	//	cout<<oldrightbase<<endl;
	//	cout<<wfmat<<endl;
	
	moveleft(leftdtmat[0],rightdtmat[sitenum-1]);
	oldleftbase=leftdtmat[0].tmatbase;
	
	//cout<<leftbase<<endl;
	//cout<<rightbase<<endl;
	//cout<<oldleftbase<<endl;
	//cout<<oldrightbase<<endl;
	//cout<<wfmat<<endl;
}

/*!
 \fn SupBlock::toComplex()
 */
void SupBlock::toComplex()
{
	if(value_type=='r')
	{
		value_type='c';
		R2C(wfmat,wfmatC);
		//cout<<wf<<endl;
		wfC.copy(wf.to_LaGenMatComplex());
		//cout<<wfC<<endl;
		for(int i=0;i<sitenum;i++)
		{
			leftdtmat[i].toComplex();
			rightdtmat[i].toComplex();
		}
	}
}
