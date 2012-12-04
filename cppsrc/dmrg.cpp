#include "dmrg.h"
#include "dtmat.h"
#include "SuperChain.h"
#include <algorithm>

namespace snake
{
namespace physics
{
    std::vector<snake::physics::Site> allfreesites;
}
}

snake::physics::DMRG::DMRG()
{
	//n=1;
	KeptStatesNum=50;
	fDMRGSweepTimes=1;
	TargetGQN.resize(1);
    std::cout<<"==========Parameters:"<<std::endl;
    std::cout<<"KeptStatesNum="<<KeptStatesNum<<std::endl;
    std::cout<<"fDMRGSweepTimes="<<fDMRGSweepTimes<<std::endl;
}


snake::physics::DMRG::~DMRG()
{
		system("rm -rf data");
}

/*!
 \fn snake::physics::DMRG::iDMRG()
 */
void snake::physics::DMRG::iDMRG()
{
	iDMRG_ReadParameters();
    std::cout<<"==========Reading free site information from matlab generated files."<<std::endl;
	iDMRG_ReadFreesites();	
    std::cout<<"==========Start Infinite DMRG."<<std::endl;
	
	///The initial left and right block each with one site.
    left=new snake::physics::Chain(snake::physics::allfreesites[0],OnSitePotentials[0]);
    right=new snake::physics::Chain(snake::physics::allfreesites[ChainLength-1],OnSitePotentials[ChainLength-1]);
	left->write("./data/L");
	right->write("./data/R");
	if(ChainLength%2==1)
	{
		///Add one more site to the initial left block if the ChainLengh is odd. 
		///This is used for the case when there is an extra impurity site at the left most.
        newleft=new snake::physics::Chain(*left,snake::physics::allfreesites[1],HoppingIntegrals[0], OnSitePotentials[1], TwoSitesInteraction[0]);
		delete left;
		left=newleft;
		left->write("./data/L");
	}
	int IniLeftL, IniRightL;
	IniLeftL=left->sitenum;
	IniRightL=right->sitenum;
	
	///Determin the occupation number for the iDMRG
    std::vector<snake::math::GQN> tempTargetGQN;
	int tempTargetGQNNum=1;
	tempTargetGQN.resize(tempTargetGQNNum);
	double occnum;
	occnum=TargetGQN[0].gqn[0]/double(ChainLength);
    std::cout<<"The average occuptation number is "<<occnum<<std::endl;
	
	///Increase both sides
	for(int i=1;i<=ChainLength/2-1;i++)
	{
		NewLeftL=IniLeftL+i;
		NewRightL=IniRightL+i;
        	tempTargetGQN[0].gqn[0]=round((NewLeftL+NewRightL)*occnum);
            //std::cout<<tempTargetGQN[0].gqn[0]<<std::endl;
		if(i>1) delete supblock;
		AddTwoSites(NewLeftL-1,NewRightL-1);
        supblock=new SuperChain(newleft,newright,left,right,H[NewLeftL-1]);
		supblock->TargetGQN=tempTargetGQN;
        //snake::math::printvector(tempTargetGQN);
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
        std::cout<<std::endl;
	}
	delete newleft;
	delete newright;
	delete supblock;
    //std::cout<<std::endl<<"The maximium exact block size is "<<n<<std::endl<<std::endl;
    std::cout<<"==========Infinite DMRG complete."<<std::endl<<std::endl;
}


/*!
 \fn snake::physics::DMRG::fDMRG()
 */
void snake::physics::DMRG::fDMRG()
{
    std::cout<<"==========Start Finite DMRG."<<std::endl;
	fDMRG_finalNewLeftL=NewLeftL;
	for(int num=0;num<fDMRGSweepTimes;num++)
	{
		//left increase at the cost of right
        std::cout<<"*****No."<<num<<" sweep"<<std::endl;
		fDMRG_Sweep2Left(NewLeftL-1,1);
		fDMRG_Sweep2Right(1,ChainLength-3);
		fDMRG_Sweep2Left(ChainLength-3,fDMRG_finalNewLeftL-1);
	}
    std::cout<<"==========Finite DMRG complete."<<std::endl<<std::endl;
}


/*!
 \fn snake::physics::DMRG::FTDMRG()
 */
void snake::physics::DMRG::FTDMRG()
{
	/// @todo implement me
}


/*!
 \fn snake::physics::DMRG::tDMRG
 */
void snake::physics::DMRG::tDMRG()
{
	
	tDMRG_ReadParameters();
    std::cout<<"==========Start t-DMRG."<<std::endl;
	supblock->creatoutputfiles();
	supblock->evolve(rt_td_impurity_OP,rt_OP,tDMRGStepNum);
	supblock->closeoutputfiles();
	delete supblock;
    std::cout<<"==========t-DMRG complete."<<std::endl<<std::endl;
}



void snake::physics::DMRG::iDMRG_ReadFreesites()
{
	std::string fname="./model/site_operators.dat";
    std::ifstream opfin(fname.c_str(),std::ios_base::in|std::ios_base::binary);
	std::string fname2="./model/site_base.dat";
    std::ifstream basefin(fname2.c_str(),std::ios_base::in|std::ios_base::binary);

    snake::physics::allfreesites.resize(ChainLength);
	for(int i=0;i<ChainLength;i++)
        snake::physics::allfreesites[i].readsite(basefin, opfin);
	opfin.close();
	basefin.close();
}

/*!
 \fn snake::physics::DMRG::ReadParameters()
 */
void snake::physics::DMRG::iDMRG_ReadParameters()
{
    std::cout<<"==========Reading model parameters from Matlab generated file."<<std::endl;
	std::string fname="./model/problemparmeters.dat";
    std::ifstream fin(fname.c_str(),std::ios_base::in|std::ios_base::binary);
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
    fin.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
	HoppingIntegrals.resize(ChainLength-1);
	OnSitePotentials.resize(ChainLength);
	TwoSitesInteraction.resize(ChainLength-1);
    snake::math::readvector(fin,HoppingIntegrals);
    snake::math::readvector(fin, OnSitePotentials);
    snake::math::readvector(fin, TwoSitesInteraction);
	fin.close();
	
	fname="./model/HC.dat";
    fin.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    snake::math::ReadOperators(fin,H,ChainLength-1);
	fin.close();

    std::cout<<"ChainLength="<<ChainLength<<std::endl;
    std::cout<<"TargetGQN="<<std::endl;
    snake::math::printvector(TargetGQN); std::cout<<std::endl;
    std::cout<<"tDMRGStepNum="<<tDMRGStepNum<<std::endl;
    std::cout<<"Hopping Integrals are:"<<std::endl;
    snake::math::printvector(HoppingIntegrals);std::cout<<std::endl;
    std::cout<<"On site potentials are:"<<std::endl;
    snake::math::printvector(OnSitePotentials);std::cout<<std::endl;
    std::cout<<"Two sites nearest neighbor interaction are:"<<std::endl;
    snake::math::printvector(TwoSitesInteraction);std::cout<<std::endl;
}


/*!
 \fn snake::physics::DMRG::mkdir
 */
void snake::physics::DMRG::mkdir()
{
	system("rm -rf results");
	system("mkdir results");
	system("rm -rf data");
	system("mkdir data");
}


/*!
 \fn snake::physics::DMRG::iDMRG2tDMRG()
 */
void snake::physics::DMRG::iDMRG2tDMRG()
{
    std::cout<<"==========Passing fDMRG supblock to the t-DMRG."<<std::endl;
	///Pass the DMRG supblock to the tDMRG
	ReadSavedLRBlocks(fDMRG_finalNewLeftL-1);
	AddTwoSites(fDMRG_finalNewLeftL-1,fDMRG_finalNewLeftL-2);
    supblock=new SuperChain(newleft,newright,left,right,H[fDMRG_finalNewLeftL-1]);
	//TargetGQN.resize(2);
	supblock->TargetGQN=TargetGQN;
	supblock->CalGroundState();
	if(newleft->base.Dim>KeptStatesNum)
		supblock->renorm(KeptStatesNum);
	supblock->KeptStatesNum=KeptStatesNum;
    //std::cout<<freesite<<std::endl;
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
    //std::cout<<wfmat1<<std::endl;std::cout<<wfmat2<<std::endl;
	for(int i=0;i<ChainLength;i++)
	{
        snake::physics::allfreesites[i].toComplex();
        snake::physics::allfreesites[i].eval();
	}
	delete left;
	delete right;
	delete newleft;
	delete newright;	
    std::cout<<std::endl;
system("rm -rf data");
    std::cout<<"==========Passing fDMRG supblock to the t-DMRG complete."<<std::endl<<std::endl;
}




/*!
 \fn snake::physics::DMRG::fDMRG_Sweep2Right(int StartChainLength, int EndChainLength)
 */
void snake::physics::DMRG::fDMRG_Sweep2Right(int StartLeftChainLength, int EndLeftChainLength)
{
    std::cout<<"-----Start sweeping to right."<<std::endl;
	for(int i=StartLeftChainLength;i<=EndLeftChainLength;i++)
	{
		ReadSavedLRBlocks(i);
		AddTwoSites(i,ChainLength-i-2);
        supblock=new SuperChain(newleft,newright,left,right,H[i]);
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
        std::cout<<std::endl;
	}
}


/*!
 \fn snake::physics::DMRG::fDMRG_Sweep2Left(int StartChainLength, int EndChainLength)
 */
void snake::physics::DMRG::fDMRG_Sweep2Left(int StartLeftChainLength, int EndLeftChainLength)
{
    std::cout<<"-----Start sweeping to left."<<std::endl;
	for(int i=StartLeftChainLength;i>=EndLeftChainLength;i--)
	{
		ReadSavedLRBlocks(i);
		AddTwoSites(i,ChainLength-i-2);
        supblock=new SuperChain(newleft,newright,left,right,H[i]);
		supblock->TargetGQN=TargetGQN;
		supblock->CalGroundState();
        //std::cout<<supblock->wf<<std::endl;break;
		//	if(newright->base.Dim>KeptStatesNum)
		supblock->renormright(KeptStatesNum);
		
		newright->write("./data/R");
		delete left;
		delete right;
		delete newright;
		delete newleft;
		delete supblock;
        std::cout<<std::endl;
	}
}


/*!
 \fn snake::physics::DMRG::ReadSavedLRBlocks(int LeftChainLength)
 */
void snake::physics::DMRG::ReadSavedLRBlocks(int LeftChainLength)
{
	std::string fname;
    std::stringstream stl,str;
	fname="./data/L";
	stl<<LeftChainLength;
	fname=fname+stl.str();
  // std::cout<<fname<<std::endl;
    left=new snake::physics::Chain(fname);
	fname="./data/R";
	str<<(ChainLength-LeftChainLength-2);
	fname=fname+str.str();
    //std::cout<<fname<<std::endl;
    right=new snake::physics::Chain(fname);
}


/*!
 \fn snake::physics::DMRG::CalN_After_fDMRG()
 */
void snake::physics::DMRG::CalN()
{
    std::cout<<"Calculate N of Everysites when Chain sweeped to the middle."<<std::endl;
	ReadSavedLRBlocks(fDMRG_finalNewLeftL-1);
	AddTwoSites(fDMRG_finalNewLeftL-1,fDMRG_finalNewLeftL-2);
    supblock=new SuperChain(newleft,newright,left,right,H[fDMRG_finalNewLeftL-1]);
	supblock->TargetGQN=TargetGQN;
	supblock->CalGroundState();
	if(newleft->base.Dim>KeptStatesNum)
		supblock->renorm(KeptStatesNum);
	
	//supblock->calN("./data/n.dat");
    std::ofstream fout("./data/n.dat");
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
 \fn snake::physics::DMRG::AddTwoSites(int LeftSitePosition)
 */
void snake::physics::DMRG::AddTwoSites(int LeftChainLength, int RightChainLength)
{
    newleft=new snake::physics::Chain(*left,snake::physics::allfreesites[LeftChainLength],HoppingIntegrals[LeftChainLength-1], OnSitePotentials[LeftChainLength], TwoSitesInteraction[LeftChainLength-1]);
    std::cout<<"NewLeftL="<<newleft->sitenum<<"\t";
    newright=new snake::physics::Chain(snake::physics::allfreesites[ChainLength-RightChainLength-1],*right,HoppingIntegrals[ChainLength-RightChainLength-1], OnSitePotentials[ChainLength-RightChainLength-1], TwoSitesInteraction[ChainLength-RightChainLength-1]);
    std::cout<<"NewRightL="<<newright->sitenum<<"\t";
}


/*!
 \fn snake::physics::DMRG::tDMRG_ReadParameters()
 */
void snake::physics::DMRG::tDMRG_ReadParameters()
{
    std::cout<<"==========Reading model parameters for t-DMRG from Matlab generated files."<<std::endl;
	std::string fname;
	fname="./model/rt_T0.dat";
    std::ifstream fin(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    snake::math::ReadOperators(fin,rt_OP,ChainLength-1);
	fin.close();
	fname="./model/rt_H1_T0.dat";
    fin.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
    snake::math::ReadOperators(fin,rt_td_impurity_OP,tDMRGStepNum);
	fin.close();
}
