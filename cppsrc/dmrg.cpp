#include "dmrg.h"
#include <algorithm>

DMRG::DMRG()
{
	//n=1;
	KeptStatesNum=50;
	fDMRGSweepTimes=1;
	TargetGQN.resize(1);
	DTMat::MaxTruncateError=1e-12;
	cout<<"==========Parameters:"<<endl;
	cout<<"KeptStatesNum="<<KeptStatesNum<<endl;
	cout<<"fDMRGSweepTimes="<<fDMRGSweepTimes<<endl;
}


DMRG::~DMRG()
{
		system("rm -rf data");
}

/*!
 \fn DMRG::iDMRG()
 */
void DMRG::iDMRG()
{
	iDMRG_ReadParameters();
	cout<<"==========Reading free site information from matlab generated files."<<endl;
	iDMRG_ReadFreesites();	
	cout<<"==========Start Infinite DMRG."<<endl;
	
	///The initial left and right block each with one site.
	left=new Block(allfreesites[0],OnSitePotentials[0]);
	right=new Block(allfreesites[ChainLength-1],OnSitePotentials[ChainLength-1]);
	left->write("./data/L");
	right->write("./data/R");
	if(ChainLength%2==1)
	{
		///Add one more site to the initial left block if the ChainLengh is odd. 
		///This is used for the case when there is an extra impurity site at the left most.
		newleft=new Block(*left,allfreesites[1],HoppingIntegrals[0], OnSitePotentials[1], TwoSitesInteraction[0]);
		delete left;
		left=newleft;
		left->write("./data/L");
	}
	int IniLeftL, IniRightL;
	IniLeftL=left->sitenum;
	IniRightL=right->sitenum;
	
	///Determin the occupation number for the iDMRG
	vector<GQN> tempTargetGQN;
	int tempTargetGQNNum=1;
	tempTargetGQN.resize(tempTargetGQNNum);
	double occnum;
	occnum=TargetGQN[0].gqn[0]/double(ChainLength);
	cout<<"The average occuptation number is "<<occnum<<endl;
	
	///Increase both sides
	for(int i=1;i<=ChainLength/2-1;i++)
	{
		NewLeftL=IniLeftL+i;
		NewRightL=IniRightL+i;
        	tempTargetGQN[0].gqn[0]=round((NewLeftL+NewRightL)*occnum);
  	        //cout<<tempTargetGQN[0].gqn[0]<<endl;
		if(i>1) delete supblock;
		AddTwoSites(NewLeftL-1,NewRightL-1);
		supblock=new SupBlock(newleft,newright,left,right,H[NewLeftL-1]);
		supblock->TargetGQN=tempTargetGQN;
		//printvector(tempTargetGQN);
		supblock->CalGroundState();
		//if(newleft->base.Dim>KeptStatesNum)
		supblock->renorm(KeptStatesNum);
		//else n++;
		newleft->write("./data/L");
		newright->write("./data/R");
		delete left;
		delete right;
		left=newleft;
		right=newright;
		cout<<endl;
	}
	delete newleft;
	delete newright;
	delete supblock;
	//cout<<endl<<"The maximium exact block size is "<<n<<endl<<endl;
	cout<<"==========Infinite DMRG complete."<<endl<<endl;
}


/*!
 \fn DMRG::fDMRG()
 */
void DMRG::fDMRG()
{
	cout<<"==========Start Finite DMRG."<<endl;
	fDMRG_finalNewLeftL=NewLeftL;
	for(int num=0;num<fDMRGSweepTimes;num++)
	{
		//left increase at the cost of right
		cout<<"*****No."<<num<<" sweep"<<endl;
		fDMRG_Sweep2Left(NewLeftL-1,1);
		fDMRG_Sweep2Right(1,ChainLength-3);
		fDMRG_Sweep2Left(ChainLength-3,fDMRG_finalNewLeftL-1);
	}
	cout<<"==========Finite DMRG complete."<<endl<<endl;
}


/*!
 \fn DMRG::FTDMRG()
 */
void DMRG::FTDMRG()
{
	/// @todo implement me
}


/*!
 \fn DMRG::tDMRG
 */
void DMRG::tDMRG()
{
	
	tDMRG_ReadParameters();
	cout<<"==========Start t-DMRG."<<endl;
	supblock->creatoutputfiles();
	supblock->evolve(rt_td_impurity_OP,rt_OP,tDMRGStepNum);
	supblock->closeoutputfiles();
	delete supblock;
	cout<<"==========t-DMRG complete."<<endl<<endl;
}



void DMRG::iDMRG_ReadFreesites()
{
	string fname="./model/site_operators.dat";
	ifstream opfin(fname.c_str(),ios_base::in|ios_base::binary);
	string fname2="./model/site_base.dat";
	ifstream basefin(fname2.c_str(),ios_base::in|ios_base::binary);

	allfreesites.resize(ChainLength);	
	for(int i=0;i<ChainLength;i++)
		allfreesites[i].readsite(basefin, opfin);
	opfin.close();
	basefin.close();
}

/*!
 \fn DMRG::ReadParameters()
 */
void DMRG::iDMRG_ReadParameters()
{
	cout<<"==========Reading model parameters from Matlab generated file."<<endl;
	string fname="./model/problemparmeters.dat";
	ifstream fin(fname.c_str(),ios_base::in|ios_base::binary);
	fin.read((char*)&ChainLength,sizeof(int));
	int TGQN;
	fin.read((char*)&TGQN,sizeof(int));
	fin.read((char*)&tDMRGStepNum,sizeof(int));
	fin.close();
	TargetGQN[0]=TGQN;
	//TargetGQN.resize(3);
	//TargetGQN[1]=TargetGQN[0].gqn[0]+1;
	//TargetGQN[2]=TargetGQN[0].gqn[0]-1;
	
	fname="./model/Hfac.dat";
	fin.open(fname.c_str(),ios_base::in|ios_base::binary);
	HoppingIntegrals.resize(ChainLength-1);
	OnSitePotentials.resize(ChainLength);
	TwoSitesInteraction.resize(ChainLength-1);
	readvector(fin,HoppingIntegrals);
	readvector(fin, OnSitePotentials);
	readvector(fin, TwoSitesInteraction);
	fin.close();
	
	fname="./model/HC.dat";
	fin.open(fname.c_str(),ios_base::in|ios_base::binary);
	ReadOperators(fin,H,ChainLength-1);
	fin.close();

	cout<<"ChainLength="<<ChainLength<<endl;
	cout<<"TargetGQN="<<endl;
	printvector(TargetGQN); cout<<endl;
	cout<<"tDMRGStepNum="<<tDMRGStepNum<<endl;
	cout<<"Hopping Integrals are:"<<endl;
	printvector(HoppingIntegrals);cout<<endl;
	cout<<"On site potentials are:"<<endl;
	printvector(OnSitePotentials);cout<<endl;
	cout<<"Two sites nearest neighbor interaction are:"<<endl;
	printvector(TwoSitesInteraction);cout<<endl;
}


/*!
 \fn DMRG::mkdir
 */
void DMRG::mkdir()
{
	system("rm -rf results");
	system("mkdir results");
	system("rm -rf data");
	system("mkdir data");
}


/*!
 \fn DMRG::iDMRG2tDMRG()
 */
void DMRG::iDMRG2tDMRG()
{
	cout<<"==========Passing fDMRG supblock to the t-DMRG."<<endl;
	///Pass the DMRG supblock to the tDMRG
	ReadSavedLRBlocks(fDMRG_finalNewLeftL-1);
	AddTwoSites(fDMRG_finalNewLeftL-1,fDMRG_finalNewLeftL-2);
	supblock=new SupBlock(newleft,newright,left,right,H[fDMRG_finalNewLeftL-1]);
	//TargetGQN.resize(2);
	supblock->TargetGQN=TargetGQN;
	supblock->CalGroundState();
	if(newleft->base.Dim>KeptStatesNum)
		supblock->renorm(KeptStatesNum);
	supblock->KeptStatesNum=KeptStatesNum;
    //cout<<freesite<<endl;
	supblock->loaddtmats();
	supblock->sweep2left(fDMRG_finalNewLeftL);
	//////////////////////////////
	supblock->TargetGQN2.resize(1);
	supblock->TargetGQN2[0].gqn[0]=supblock->TargetGQN[0].gqn[0];
	
	
	///Generate new wave functions
	Rmatrix wfmat2=supblock->wfmat;
	int row=wfmat2.rowbase.subnum;
	int col=wfmat2.colbase.subnum;
	for(int i=0;i<row;i++)
		if(wfmat2.pmat(i,0)!=-1)
		{
			wfmat2.pmat(i,1)=wfmat2.pmat(i,0);
			wfmat2.pmat(i,0)=-1;
		}
    supblock->extractwf(supblock->wfmat,supblock->wf,supblock->TargetGQN);
	supblock->toComplex();
	R2C(wfmat2,supblock->wfmatC2);  supblock->extractwf(supblock->wfmatC2,supblock->wfC2,supblock->TargetGQN2);
	//cout<<wfmat1<<endl;cout<<wfmat2<<endl;
	for(int i=0;i<ChainLength;i++)
	{
		allfreesites[i].toComplex();
		allfreesites[i].eval();
	}
	delete left;
	delete right;
	delete newleft;
	delete newright;	
	cout<<endl;
system("rm -rf data");
	cout<<"==========Passing fDMRG supblock to the t-DMRG complete."<<endl<<endl;
}




/*!
 \fn DMRG::fDMRG_Sweep2Right(int StartChainLength, int EndChainLength)
 */
void DMRG::fDMRG_Sweep2Right(int StartLeftChainLength, int EndLeftChainLength)
{
	cout<<"-----Start sweeping to right."<<endl;
	for(int i=StartLeftChainLength;i<=EndLeftChainLength;i++)
	{
		ReadSavedLRBlocks(i);
		AddTwoSites(i,ChainLength-i-2);
		supblock=new SupBlock(newleft,newright,left,right,H[i]);
		supblock->TargetGQN=TargetGQN;
		supblock->CalGroundState();
		//	if(newleft->base.Dim>KeptStatesNum)
		supblock->renormleft(KeptStatesNum);
		newleft->write("./data/L");
		delete left;
		delete right;
		delete newleft;
		delete newright;
		delete supblock;
		cout<<endl;
	}
}


/*!
 \fn DMRG::fDMRG_Sweep2Left(int StartChainLength, int EndChainLength)
 */
void DMRG::fDMRG_Sweep2Left(int StartLeftChainLength, int EndLeftChainLength)
{
	cout<<"-----Start sweeping to left."<<endl;
	for(int i=StartLeftChainLength;i>=EndLeftChainLength;i--)
	{
		ReadSavedLRBlocks(i);
		AddTwoSites(i,ChainLength-i-2);
		supblock=new SupBlock(newleft,newright,left,right,H[i]);
		supblock->TargetGQN=TargetGQN;
		supblock->CalGroundState();
		//cout<<supblock->wf<<endl;break;
		//	if(newright->base.Dim>KeptStatesNum)
		supblock->renormright(KeptStatesNum);
		
		newright->write("./data/R");
		delete left;
		delete right;
		delete newright;
		delete newleft;
		delete supblock;
		cout<<endl;
	}
}


/*!
 \fn DMRG::ReadSavedLRBlocks(int LeftChainLength)
 */
void DMRG::ReadSavedLRBlocks(int LeftChainLength)
{
	string fname;
	stringstream stl,str;
	fname="./data/L";
	stl<<LeftChainLength;
	fname=fname+stl.str();
  // cout<<fname<<endl;
	left=new Block(fname);
	fname="./data/R";
	str<<(ChainLength-LeftChainLength-2);
	fname=fname+str.str();
	//cout<<fname<<endl;
	right=new Block(fname);
}


/*!
 \fn DMRG::CalN_After_fDMRG()
 */
void DMRG::CalN()
{
	cout<<"Calculate N of Everysites when Chain sweeped to the middle."<<endl;
	ReadSavedLRBlocks(fDMRG_finalNewLeftL-1);
	AddTwoSites(fDMRG_finalNewLeftL-1,fDMRG_finalNewLeftL-2);
	supblock=new SupBlock(newleft,newright,left,right,H[fDMRG_finalNewLeftL-1]);
	supblock->TargetGQN=TargetGQN;
	supblock->CalGroundState();
	if(newleft->base.Dim>KeptStatesNum)
		supblock->renorm(KeptStatesNum);
	
	//supblock->calN("./data/n.dat");
	ofstream fout("./data/n.dat");
	newleft->calN(fout,'l');
	newright->calN(fout,'r');
	fout.close();
	
	system("xmgrace ./data/n.dat &");
	//supblock->prepare(); ///Store base information of the blocks
	delete left;
	delete right;
	delete newleft;
	delete newright;
	delete supblock;
}

/*!
 \fn DMRG::AddTwoSites(int LeftSitePosition)
 */
void DMRG::AddTwoSites(int LeftChainLength, int RightChainLength)
{
	newleft=new Block(*left,allfreesites[LeftChainLength],HoppingIntegrals[LeftChainLength-1], OnSitePotentials[LeftChainLength], TwoSitesInteraction[LeftChainLength-1]);
	cout<<"NewLeftL="<<newleft->sitenum<<"\t";
	newright=new Block(allfreesites[ChainLength-RightChainLength-1],*right,HoppingIntegrals[ChainLength-RightChainLength-1], OnSitePotentials[ChainLength-RightChainLength-1], TwoSitesInteraction[ChainLength-RightChainLength-1]);
    cout<<"NewRightL="<<newright->sitenum<<"\t";
}


/*!
 \fn DMRG::tDMRG_ReadParameters()
 */
void DMRG::tDMRG_ReadParameters()
{
	cout<<"==========Reading model parameters for t-DMRG from Matlab generated files."<<endl;
	string fname;
	fname="./model/rt_T0.dat";
	ifstream fin(fname.c_str(),ios_base::in|ios_base::binary);
	ReadOperators(fin,rt_OP,ChainLength-1);
	fin.close();
	fname="./model/rt_H1_T0.dat";
	fin.open(fname.c_str(),ios_base::in|ios_base::binary);
	ReadOperators(fin,rt_td_impurity_OP,tDMRGStepNum);
	fin.close();
}
