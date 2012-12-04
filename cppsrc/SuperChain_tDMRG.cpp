#include "SuperChain.h"


void snake::physics::SuperChain::evolve(std::vector<LaGenMatComplex> &rt_td_impurity_OP, std::vector<LaGenMatComplex> &rt_OP, int timesteps)
{
	/////////////////////Output reduced density matrix////////////////////////
	std::string filename;
	filename="./results/rdm.dat";
	std::ofstream fout_rdm(filename.c_str());
	fout_rdm.width(20);
	fout_rdm.precision(15);
	//////////////////////////////////////////////////////////////////////////
	
	COMPLEX ndot;
	
  //genindex();
	//std::cout<<wfC<<std::endl;
	//std::cout<<wfmatC<<std::endl;
	normalize(wfC);
	normalize(wfC2);
	//std::cout<<wfC<<std::endl;
	
	///Time evolve
	for(int n=0;n<timesteps;n++)
	{
		std::cout<<"TimeStep="<<n<<std::endl;
		onetimestep(rt_td_impurity_OP[n],rt_OP);
		///vonNeumann entropy when there are two sites on the left
		for(int i=2;i<sitenum-2;i++)
		{
			fout_entropy_t<<rightdtmat[sitenum-i].entropy<<" ";
		}
		fout_entropy_t<<std::endl;
		
		evalwfmat(wfC,wfmatC,TargetGQN);
#if TARGET_TWO_WF
		evalwfmat(wfC2,wfmatC2,TargetGQN2);
#endif
		//std::cout<<wfmatC<<std::endl; std::cout<<wfmatC2<<std::endl;
		
		Cmatrix ReducedDM(leftbase,leftbase);
		Mat_Trans_Mat_Mult(wfmatC, wfmatC, ReducedDM);
		//std::cout<<ReducedDM<<std::endl;
#if TARGET_TWO_WF
		
		Cmatrix ReducedDM2(leftbase,leftbase);
		Mat_Trans_Mat_Mult(wfmatC2, wfmatC2, ReducedDM2);
		ReducedDM+=ReducedDM2;
		//std::cout<<ReducedDM2<<std::endl;
		//Cmatrix ReducedDM3(leftbase,leftbase);
		//std::cout<<wfmatC2<<std::endl;
		//std::cout<<wfmatC<<std::endl;
		//std::cout<<ReducedDM3<<std::endl;
		ReducedDM2=0;
		Mat_Trans_Mat_Mult(wfmatC2, wfmatC, ReducedDM2);
		//std::cout<<ReducedDM2<<std::endl;
		ReducedDM+=ReducedDM2;
		
		ReducedDM2=0;
		Mat_Trans_Mat_Mult(wfmatC, wfmatC2, ReducedDM2);
		ReducedDM+=ReducedDM2;
		//std::cout<<ReducedDM<<std::endl;
		//fout_rdm<<ReducedDM.submat[ReducedDM.pmat(0,0)](0,0)<<" ";
		//fout_rdm<<ReducedDM.submat[ReducedDM.pmat(1,0)](0,0)<<" ";
		//fout_rdm<<ReducedDM.submat[ReducedDM.pmat(0,1)](0,0)<<" ";
		//fout_rdm<<ReducedDM.submat[ReducedDM.pmat(1,1)](0,0)<<std::endl;
		
		fout_rdm<<ReducedDM.submat[ReducedDM.pmat(1,0)](0,0).r<<" ";
		fout_rdm<<ReducedDM.submat[ReducedDM.pmat(1,0)](0,0).i<<std::endl;
		
		
#endif
  	//std::cout<<ReducedDM<<std::endl;
		ndot=ReducedDM.submat[ReducedDM.pmat(0,0)](0,0);
	 	COMPLEX offdiag1=ReducedDM.submat[ReducedDM.pmat(0,0)](0,1);
		COMPLEX offdiag2=ReducedDM.submat[ReducedDM.pmat(0,0)](1,0);
		//std::cout<<"offdiag1="<<offdiag1<<"\t|"<<"offdiag2="<<offdiag2<<"\t";
		double sigmaz=1-2*ndot.r;
		fout_1stsite_n_t<<sigmaz<<std::endl;
		std::cout<<"sigma_z(t)="<<sigmaz<<"\t";
		fout_rdm<<offdiag1.r<<" "<<offdiag1.i<<std::endl;
		std::cout<<std::endl<<std::endl;
		
	}
  //delindex();
}

void snake::physics::SuperChain::onetimestep(LaGenMatComplex& impurity_OP_t, std::vector<LaGenMatComplex>& rt_OP)
{
	///Move right
	//std::cout<<"LKeptStatesNum="<<flush;
	for(int i=1;i<=sitenum-1;i++)
	{
		if(i==1)
		{
  		//std::cout<<impurity_OP_t<<std::endl;
  		//std::cout<<TargetGQN[0]<<std::endl;
        midsite1=snake::physics::allfreesites[0];
        midsite2=snake::physics::allfreesites[1];
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
        midsite1=snake::physics::allfreesites[i-1];
        midsite2=snake::physics::allfreesites[i];
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
            freesite=snake::physics::allfreesites[i+1];
			moveright(leftdtmat[i],rightdtmat[sitenum-i-1]);
			oldrightbase=rightdtmat[sitenum-i-2].tmatbase;
		}
		
	}
	std::cout<<std::endl;
	
	
	///Move left
	//std::cout<<"RKeptStatesNum="<<flush;
	for(int i=sitenum-1;i>=1;i--)
	{
		if(i==1)
		{
  		//std::cout<<TargetGQN[0]<<std::endl;
        midsite1=snake::physics::allfreesites[0];
        midsite2=snake::physics::allfreesites[1];
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
        midsite1=snake::physics::allfreesites[i-1];
        midsite2=snake::physics::allfreesites[i];
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
            freesite=snake::physics::allfreesites[i];
			moveleft(leftdtmat[i-1],rightdtmat[sitenum-i]);
			oldleftbase=leftdtmat[i-2].tmatbase;
		}
		
	}
	std::cout<<std::endl;
	normalize(wfC);
	normalize(wfC2);
}


void snake::physics::SuperChain::creatoutputfiles()
{
	std::string filename;
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

void snake::physics::SuperChain::closeoutputfiles()
{
#if CAL_DURING_TIME_EVOLVE
	fout_1stsite_n_t.close();
#endif
	fout_entropy_t.close();
	fout_steperror_t.close();
	
}


void snake::physics::SuperChain::loaddtmats()
{
	rightdtmat.resize(sitenum);
	leftdtmat.resize(sitenum);
	
	
	for(int i=1;i<=sitenum/2;i++)
	{
		std::string filename;
        std::stringstream stemp2;
		stemp2<<i;
		filename="./data/L";
		filename+=stemp2.str();
        snake::physics::Chain lb(filename);
		filename="./data/R";
		filename+=stemp2.str();
        snake::physics::Chain rb(filename);
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
	
	
	//std::cout<<leftdtmat[1].tmatbase<<std::endl;
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


void snake::physics::SuperChain::sweep2left(int NewLeftChainLength)
{
	for(int i=NewLeftChainLength;i>1;i--)
	{
		//std::cout<<"Left site is "<<i<<std::endl;
		//std::cout<<wfmat<<std::endl;
		//std::cout<<"Norm="<<Blas_Norm2(wf)<<std::endl;
            freesite=snake::physics::allfreesites[i];
		moveleft(leftdtmat[i-1],rightdtmat[sitenum-i]);
		oldleftbase=leftdtmat[i-2].tmatbase;
		
	}
}

void snake::physics::SuperChain::sweep2leftmost()
{
	
	//	std::cout<<leftbase<<std::endl;
	//	std::cout<<rightbase<<std::endl;
	//	std::cout<<oldleftbase<<std::endl;
	//	std::cout<<oldrightbase<<std::endl;
	//	std::cout<<wfmat<<std::endl;
	
	moveleft(leftdtmat[0],rightdtmat[sitenum-1]);
	oldleftbase=leftdtmat[0].tmatbase;
	
	//std::cout<<leftbase<<std::endl;
	//std::cout<<rightbase<<std::endl;
	//std::cout<<oldleftbase<<std::endl;
	//std::cout<<oldrightbase<<std::endl;
	//std::cout<<wfmat<<std::endl;
}


void snake::physics::SuperChain::toComplex()
{
	if(value_type=='r')
	{
		value_type='c';
		R2C(wfmat,wfmatC);
		//std::cout<<wf<<std::endl;
		wfC.copy(wf.to_LaGenMatComplex());
		//std::cout<<wfC<<std::endl;
		for(int i=0;i<sitenum;i++)
		{
			leftdtmat[i].toComplex();
			rightdtmat[i].toComplex();
		}
	}
}
