/*
 *  supblock_FTDMRG.cpp
 *  snake
 *
 *  Created by Cheng Guo on 12/14/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include"SuperChain.h"

/*!
 \fn snake::physics::SuperChain::geninitialwf();
 */
void snake::physics::SuperChain::geninitialwf()
{
	wfmat.resize(rightbase,leftbase);
	wfmat.subnum=1;
	wfmat.pmat(1,1)=wfmat.subnum-1;
	wfmat.submat.resize(wfmat.subnum);
	wfmat.submat[wfmat.subnum-1].resize(2,2);
	wfmat.submat[wfmat.subnum-1]=0;
	wfmat.submat[wfmat.subnum-1](0,0)=1;
}


/*!
 \fn snake::physics::SuperChain::geninitialbases()
 */
void snake::physics::SuperChain::geninitialbases()
{
	leftbase.Dim=4;
	leftbase.subnum=3;
	leftbase.dim.resize(leftbase.subnum);
	leftbase.subgqn.resize(leftbase.subnum);
	leftbase.dim[0]=1;
	leftbase.dim[1]=2;
	leftbase.dim[2]=1;
	leftbase.subgqn[0]=-1;
	leftbase.subgqn[1]=0;
	leftbase.subgqn[2]=1;

	rightbase=leftbase;
	oldrightbase.Dim=1;
	oldrightbase.subnum=1;
	oldrightbase.dim.resize(oldrightbase.subnum);
	oldrightbase.subgqn.resize(oldrightbase.subnum);
	oldrightbase.dim[0]=1;
	oldrightbase.subgqn[0]=0;

	oldleftbase=oldrightbase;
}


/*!
 \fn snake::physics::SuperChain::geninitialdtmats(std::vector<DTMat> &left, std::vector<DTMat> &right)
*/
void snake::physics::SuperChain::geninitialdtmats(std::vector<DTMat> &leftdt, std::vector<DTMat> &rightdt)
{
	LaGenMatDouble trunmat,trunmattrans;
	leftdt.resize(sitenum);
	rightdt.resize(sitenum);
	for(int i=0;i<sitenum;i++)
	{
		leftdt[i].handside='l';
		rightdt[i].handside='r';
		rightdt[i].tmatbase.Dim=1;
		rightdt[i].tmatbase.subnum=1;
		rightdt[i].tmatbase.dim.resize(1);
		rightdt[i].tmatbase.subgqn.resize(1);
		rightdt[i].tmatbase.dim[0]=1;
		rightdt[i].tmatbase.subgqn[0]=0;
		rightdt[i].rightbase=freesite.base;

		trunmat.resize(freesite.base.Dim,1);
		trunmattrans.resize(1,freesite.base.Dim);
		trunmattrans=0;
		trunmat=0;
		trunmat(1,0)=1;
		trunmattrans(0,1)=1;
		rightdt[i].trunmat=Rmatrix(trunmat,freesite.base,rightdt[i].tmatbase,1);
		//std::cout<<rightdt[i].trunmat;
		//rightdt[i].trunmattrans=Rmatrix(trunmattrans,rightdt[i].tmatbase,freesite.base,1);
		// std::cout<<rightdt[i].trunmattrans;
	}
	leftdt[0].tmatbase=rightdt[0].tmatbase;

}

/*!
 \fn snake::physics::SuperChain::evolve(LaGenMatDouble &T,int timesteps)
 */
void snake::physics::SuperChain::evolve(std::vector<LaGenMatDouble> &gs_proj,int timesteps)
{
	///Initialize
	geninitialdtmats(leftdtmat,rightdtmat);
	geninitialbases();
	geninitialwf();
	genindex();

	LaGenMatDouble H;
	LaVectorDouble tempwf;
	GQNBase temprightbase;

	LaVectorDouble energy(timesteps);
	energy=0;

	///Time evolve
	extractwf(wfmat,wf, TargetGQN);
	for(int n=0;n<timesteps;n++)
	{
		std::cout<<std::endl<<"********************Temperature step "<<n<<"********************"<<std::endl;

		///Move right
		for(int i=1;i<sitenum-1;i++)
		{
			//std::cout<<"Left site is "<<i<<std::endl;

			genmiddlemap(TargetGQN);
			middlemult(gs_proj[i-1],wf);
			deletemiddlemap();
			moveright(leftdtmat[i],rightdtmat[sitenum-i-1]);
			oldrightbase=rightdtmat[sitenum-i-2].tmatbase;
		}
		genmiddlemap(TargetGQN);
		middlemult(gs_proj[sitenum-2],wf);
		deletemiddlemap();
		//wf.scale(1/sqrt(Blas_Dot_Prod(wf,wf)));
		//chop(wf,1e-15);std::cout<<wf<<std::endl;

		///Move left
		for(int i=sitenum-1;i>1;i--)
		{
			//std::cout<<"Left site is "<<i<<std::endl;
			genmiddlemap(TargetGQN);
			middlemult(gs_proj[i-1],wf);
			deletemiddlemap();
			moveleft(leftdtmat[i-1],rightdtmat[sitenum-i]);
			if(i==10)
			{
				//fentropy<<rightdtmat[sitenum-i].entropy<<std::endl;
				std::cout<<"The von Neumann entropy of this step is "<<rightdtmat[sitenum-i].entropy<<std::endl;
			}
			oldleftbase=leftdtmat[i-2].tmatbase;
		}
		genmiddlemap(TargetGQN);
		middlemult(gs_proj[0],wf);
		deletemiddlemap();
		double mod=sqrt(Blas_Dot_Prod(wf,wf));
		std::cout<<"The mod of wf is "<<mod<<std::endl;
		wf.scale(1/sqrt(Blas_Dot_Prod(wf,wf)));


		//Caculate on site occupation number
		evalwfmat(wf,wfmat,TargetGQN);
		Rmatrix leftdm(leftbase,leftbase);
		Mat_Trans_Mat_Mult(wfmat,wfmat,leftdm);
		double occnum;
		Rmatrix tempmat(leftbase,leftbase);
		Mat_Mat_Mult(leftdm,freesite.n[0],tempmat);
		//std::cout<<leftdm<<std::endl;
		//std::cout<<freesite.n[0]<<std::endl;
		occnum=tempmat.trace();
		//    foccnum<<occnum<<std::endl;
	}
	//foccnum<<std::endl;
	//fentropy<<std::endl;
	delindex();
	//system("xmgrace ./data/entropy_temperature_evolve &");
}
