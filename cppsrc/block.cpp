#include "block.h"
#include "dtmat.h"
#include "blocham.h"
#include "site.h"
#include "public.h"
#include "setting.h"

#include <blaspp.h>
#include <sstream>
#include <iostream>
using namespace std;

Block::Block()
{
	value_type='r';
}

Block::~Block()
{
	delete hamiltonian;

	delete dtmat;
}

/*!
 \fn Block::Block(Site site)
 */
Block::Block(Site &first, double OnSiteE)
{
	value_type=first.value_type;
	hamiltonian=new BlocHam(first, OnSiteE);
	site.resize(1);
	site[0]=first;
	siteadded=&first;
	sitenum=1;
	dtmat=new DTMat;
	base=first.base;
}

/*!
 \fn Block::Block(char *filename)
 */
Block::Block(std::string &filename)
{
	read(filename);
}


/*!
 \fn Block::Block(Block &old,Site add)
 */
/*
 Block::Block(Block &old,Site &addsite)
 {
 value_type=old.value_type;
 #if FERMIONSIGN
 if(old.sitenum==1) old.site[0].multsignmat();
 #endif
 if(old.sitenum==1)
 {
 old.hamiltonian->H=old.site[0].n[0];
 old.hamiltonian->H.scale(-10000.0);
 std::cout<<"!!!!!!!!!!!!!!!!Be Careful about the potential on the impurity!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
 std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Edit block.cpp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
 std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!This is a alpha version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
 }
 //old.site[0]
 initialadd(old,addsite);
 calHinter(old.site[old.sitenum-1],addsite,'r');
 if(value_type=='r')
 {
 hamiltonian=new BlocHam(old.hamiltonian,addsite,Hinter,'r');
 }
 else
 {
 hamiltonian=new BlocHam(old.hamiltonian,addsite,HinterC,'r');
 }

 base=hamiltonian->base;

 #if CALCULATE_ALL_SITES
 ///Sites in old block
 for(int i=0;i<old.sitenum;i++)
 {
 site[i]=old.site[i];
 site[i].addsite(addsite,'r');
 }
 #endif

 ///The site just added
 site[sitenum-1]=addsite;
 site[sitenum-1].addtoblock(old,'r');

 dtmat=new DTMat(value_type,'l');
 std::cout<<"Block dim is "<<base.Dim<<std::endl;
 }

 Block::Block(Site &addsite,Block &old)
 {
 value_type=old.value_type;

 initialadd(old,addsite);
 calHinter(addsite,old.site[old.sitenum-1],'l');
 if(value_type=='r')
 {
 hamiltonian=new BlocHam(old.hamiltonian,addsite,Hinter,'l');
 }
 else
 {
 hamiltonian=new BlocHam(old.hamiltonian,addsite,HinterC,'l');

 }

 base=hamiltonian->base;
 #if CALCULATE_ALL_SITES
 for(int i=0;i<old.sitenum;i++)
 {
 site[i]=old.site[i];
 site[i].addsite(addsite,'l');
 }
 #endif
 site[old.sitenum]=addsite;
 site[old.sitenum].addtoblock(old,'l');
 dtmat=new DTMat(value_type,'r');
 std::cout<<"Block dim is "<<base.Dim<<std::endl;
 }
 */

Block::Block(Block &old,Site &addsite, double HoppingT, double OnSiteE, double TwoSitesV)
{
	value_type=old.value_type;
#if FERMIONSIGN
	if(old.sitenum==1) old.site[0].multsignmat();
#endif
	initialadd(old,addsite);
	calHinter(old.site[old.sitenum-1],addsite,'r',HoppingT, OnSiteE,TwoSitesV);
	if(value_type=='r')
	{
     		hamiltonian=new BlocHam(old.hamiltonian,addsite,Hinter,'r');
	}
	else
	{
		hamiltonian=new BlocHam(old.hamiltonian,addsite,HinterC,'r');
	}
	base=hamiltonian->base;

#if CALCULATE_ALL_SITES
	///Sites in old block
	for(int i=0;i<old.sitenum;i++)
	{
		site[i]=old.site[i];
		site[i].addsite(addsite,'r');
	}
#endif

	///The site just added
	site[sitenum-1]=addsite;
	site[sitenum-1].addtoblock(old,'r');

	dtmat=new DTMat(value_type,'l');
}


Block::Block(Site &addsite,Block &old,double HoppingT, double OnSiteE, double TwoSitesV)
{
	value_type=old.value_type;

	initialadd(old,addsite);
	calHinter(addsite,old.site[old.sitenum-1],'l',HoppingT, OnSiteE,TwoSitesV);
	if(value_type=='r')
	{
		hamiltonian=new BlocHam(old.hamiltonian,addsite,Hinter,'l');
	}
	else
	{
		hamiltonian=new BlocHam(old.hamiltonian,addsite,HinterC,'l');

	}

	base=hamiltonian->base;
#if CALCULATE_ALL_SITES
	for(int i=0;i<old.sitenum;i++)
	{
		site[i]=old.site[i];
		site[i].addsite(addsite,'l');
	}
#endif
	site[old.sitenum]=addsite;
	site[old.sitenum].addtoblock(old,'l');
	dtmat=new DTMat(value_type,'r');
}

/*
 Block::Block(Site &addsite,Block &old)
 {
 value_type=old.value_type;

 initialadd(old,addsite);
 calHinter(addsite,old.site[old.sitenum-1],'l');
 if(value_type=='r')
 {
 hamiltonian=new BlocHam(old.hamiltonian,addsite,Hinter,'l');
 }
 else
 {
 hamiltonian=new BlocHam(old.hamiltonian,addsite,HinterC,'l');

 }

 base=hamiltonian->base;
 #if CALCULATE_ALL_SITES
 for(int i=0;i<old.sitenum;i++)
 {
 site[i]=old.site[i];
 site[i].addsite(addsite,'l');
 }
 #endif
 site[old.sitenum]=addsite;
 site[old.sitenum].addtoblock(old,'l');
 dtmat=new DTMat(value_type,'r');
 std::cout<<"Block dim is "<<base.Dim<<std::endl;
 }
 */

/*!
 \fn Block::Block(Site& addsite,Block& old,int localsite)
 */
/*
Block::Block(Site& addsite,Block& old,int localsite)
{
	value_type=old.value_type;
	initialadd(old,addsite);

	if(value_type=='r')
	{
		if(localsite==sitenum-1)
		{
			calHinter(addsite,old.site[old.sitenum-1],'l');
			hamiltonian=new BlocHam(old.hamiltonian,addsite,Hinter,'l');
		}
		else
			hamiltonian=new BlocHam(old.hamiltonian,addsite,'l');
	}
	else
	{
		if(localsite==sitenum-1)
		{
			calHinter(addsite,old.site[old.sitenum-1],'l');
			hamiltonian=new BlocHam(old.hamiltonian,addsite,HinterC,'l');
		}
		else
			hamiltonian=new BlocHam(old.hamiltonian,addsite,'l');
	}


	base=hamiltonian->base;
#if CALCULATE_ALL_SITES
	for(int i=0;i<old.sitenum;i++)
	{
		site[i]=old.site[i];
		site[i].addsite(addsite,'l');
	}
#endif
	site[old.sitenum]=addsite;
	site[old.sitenum].addtoblock(old,'l');


	dtmat=new DTMat(value_type,'r');
}
 */

/*!
 \fn Block::initialadd(Block &old,Site &add)
 */
void Block::initialadd(Block &old,Site &add)
{
	siteadded=&add;///Might be rewrite the way to store steadded.
	sitenum=old.sitenum+1;
	site.resize(sitenum);
}


/*!
 \fn Block::interaction(Site &a,Site &b)
 */
/*
 void Block::calHinter(Site &siteA,Site &siteB,char addsiteposition,char include_onsite)
 {
 calHinter_Heisenberg(siteA,siteB,addsiteposition);
 //calHinter_Hubbard(siteA,siteB,addsiteposition,include_onsite);
 //std::cout<<"Hinter"<<Hinter<<std::endl;

 }
 */
void Block::calHinter(Site &siteA,Site &siteB,char addsiteposition,double HoppingT, double OnSiteE, double TwoSitesV)
{
	//std::cout<<"HoppingT="<<HoppingT<<std::endl;


	if(value_type=='r')
	{
		//Hopping term
  	GQNBase Hinterbase;
  	Hinterbase=kron(siteA.base,siteB.base);
  	
	 Hinter.geneye(Hinterbase);
	Hinter=0;
	 if(HoppingT!=0)
	 {  	 
  	        Rmatrix tempx,tempy;
		siteA.eval();
		siteB.eval();

		kron(siteA.c[0],siteB.a[0],tempx);
		kron(siteA.a[0],siteB.c[0],tempy);
		tempx+=tempy;
		tempx.scale(HoppingT);
  	 
  	       Hinter+=tempx;
  	}
		//Onsite Potential Part
		if(OnSiteE!=0)
		{
			std::cout<<"OnsiteE="<<OnSiteE<<std::endl;
			Rmatrix tempeye,H_onsite;
			if(addsiteposition=='r')
			{
				tempeye.geneye(siteA.base);
				kron(tempeye,siteB.n[0],H_onsite);
			}
			else {
				tempeye.geneye(siteB.base);
				kron(siteA.n[0],tempeye,H_onsite);
			}
			H_onsite.scale(OnSiteE);
			Hinter+=H_onsite;
		}
  	//std::cout<<Hinter<<std::endl;

		//Two sites interaction term
		//std::cout<<"!!!!!!!!!!Model Dependent Part!!!!!!!!!!"<<std::endl;
		//std::cout<<"H_i=V(n_i-1/2)(n_{i+1}-1/2)"<<std::endl;
		if(TwoSitesV!=0)
		{
			std::cout<<"TwoSitesV="<<TwoSitesV<<std::endl;
			Rmatrix temp1,temp2,temp3;
			temp1.geneye(siteA.base);
			temp2.geneye(siteB.base);
			temp1.scale(-0.5);
			temp2.scale(-0.5);
			temp1=siteA.n[0]+temp1;
			temp2=siteB.n[0]+temp2;
			kron(temp1,temp2,temp3);
			temp3.scale(TwoSitesV);
			Hinter+=temp3;
		}

	}
	else
	{
		std::cout<<"calHinter_Heisenberg() has not been inplemented for complex numbers"<<std::endl;
	}
  //std::cout<<Hinter<<std::endl;
}

/*!
 \fn Block::renorm(LaGenMatDouble &tmat)
 */
void Block::renorm()
{
	hamiltonian->renorm(*dtmat);

#if CALCULATE_ALL_SITES
	renormsites();
#else
	renormsidesite();
#endif
	base=dtmat->tmatbase;
}



/*!
 \fn Block::renormsites(DTMat &mat)
 */
void Block::renormsites()
{
	for(int i=0;i<sitenum;i++)
		site[i].renorm(*dtmat);
}


/*!
 \fn Block::renormsidesite(DTMat &mat)
 */
void Block::renormsidesite()
{
	site[sitenum-1].renorm(*dtmat);
}


/*!
 \fn Block::write(char *prefix)
 */
void Block::write(char *prefix)
{
	std::string filename;
	stringstream stemp;
	stemp<<sitenum;
	filename=prefix+stemp.str();
  //std::cout<<std::endl<<filename<<std::endl;
	std::ofstream fout(filename.c_str(),std::ios_base::out|std::ios_base::binary);
	fout.write(&value_type,sizeof value_type);
	fout.write((char*)&sitenum,sizeof sitenum);
#if CALCULATE_ALL_SITES
	for(int i=0;i<sitenum;i++)    site[i].write(fout);
#else
	site[sitenum-1].write(fout);
#endif
	hamiltonian->write(fout);

	dtmat->write(fout);
	base.write(fout);
	fout.close();
}


/*!
 \fn Block::read(char *filename)
 */
void Block::read(std::string &filename)
{
	std::ifstream fin(filename.c_str(),std::ios_base::in|std::ios_base::binary);
	fin.read(&value_type,sizeof value_type);
	fin.read((char*)&sitenum,sizeof sitenum);
	site.resize(sitenum);
#if CALCULATE_ALL_SITES
	for(int i=0;i<sitenum;i++) site[i].read(fin);
#else
	site[sitenum-1].read(fin);
#endif
	hamiltonian=new BlocHam(fin);

	dtmat=new DTMat(fin);
	base.read(fin);
	fin.close();
}


/*!
 \fn Block::calN(std::ofstream fout)
 */
void Block::calN(std::ofstream &fout,char hand)
{
	if(value_type=='r')
	{
		Rmatrix tempmat(base,base);
		dtmat->renorm();
		if(hand=='l')
			for(int i=0;i<sitenum;i++)
			{
				//std::cout<<"Calculate Site "<<i<<std::endl;
				//std::cout<<dtmat->denmat.size(0)<<std::endl;
				//std::cout<<site[i].n.size(0)<<std::endl;
				for(int j=0;j<site[i].num;j++)
				{
					Mat_Mat_Mult(dtmat->denmat,site[i].n[j],tempmat);
					fout<<tempmat.trace()<<" ";
				}
				fout<<std::endl;
				//std::cout<<trace(tempmat)<<std::endl;
			}
		else
			for(int i=sitenum-1;i>=0;i--)
			{
				for(int j=0;j<site[i].num;j++)
				{
					Mat_Mat_Mult(dtmat->denmat,site[i].n[j],tempmat);
					fout<<tempmat.trace()<<" ";
				}
				fout<<std::endl;
			}
	}
	else
	{
		Cmatrix tempmat2(base,base);
		Rmatrix tempmat(base,base);
		dtmat->renorm();
		if(hand=='l')
			for(int i=0;i<sitenum;i++)
			{
				//std::cout<<"Calculate Site "<<i<<std::endl;
				//std::cout<<dtmat->denmat.size(0)<<std::endl;
				//std::cout<<site[i].n.size(0)<<std::endl;
				for(int j=0;j<site[i].num;j++)
				{
					Mat_Mat_Mult(dtmat->denmatC,site[i].nC[j],tempmat2);
					C2R(tempmat2,tempmat);
					fout<<tempmat.trace()<<" ";
				}
				std::cout<<std::endl;
				//std::cout<<trace(tempmat)<<std::endl;
			}
		else
			for(int i=sitenum-1;i>=0;i--)
			{
				for(int j=0;j<site[i].num;j++)
				{
					Mat_Mat_Mult(dtmat->denmatC,site[i].nC[j],tempmat2);
					C2R(tempmat2,tempmat);
					fout<<tempmat.trace()<<" ";
				}
				std::cout<<std::endl;
			}
	}

}


/*!
 \fn Block::toComplex()
 */
void Block::toComplex()
{
	if(value_type=='r')
	{
		value_type='c';
		hamiltonian->toComplex();
		dtmat->toComplex();
		for(int i=0;i<sitenum;i++)
			site[i].toComplex();
		siteadded->toComplex();
	}
}


/*!
 \fn Block::calHinter_Hubbard(Site &siteA,Site &siteB)
 */
/*
void Block::calHinter_Hubbard(Site &siteA,Site &siteB,char addsiteposition,char include_onsite)
{

	if(value_type=='r')
	{
		Rmatrix temp,temp1,temp2,tempeye;
		siteA.eval();
		siteB.eval();

		if(addsiteposition=='r')
		{
			//Hopping term
			//std::cout<<siteA.a[0]<<std::endl;
			//std::cout<<siteB.c[0]<<std::endl;
			kron(siteA.a[0],siteB.c[0],temp1);
			//std::cout<<temp1<<std::endl;
			kron(siteA.a[1],siteB.c[1],temp);
			//std::cout<<temp<<std::endl;
			temp+=temp1;
			//std::cout<<temp<<std::endl;
			//std::cout<<temp.rowbase<<std::endl;
			//std::cout<<temp.colbase<<std::endl;
			Transpose(temp,temp2);
			//std::cout<<temp<<std::endl;
			//std::cout<<temp2<<std::endl;
			temp+=temp2;
			temp.scale(T);
			//std::cout<<temp<<std::endl;
		}
		else
		{
			//Hopping term
			Rmatrix signmat;
			signmat.gensignmat(siteA.base,0);
			//std::cout<<signmat<<std::endl;
			temp=siteA.a[0]*signmat;
			kron(temp,siteB.c[0],temp1);
			signmat.gensignmat(siteA.base,1);
			//std::cout<<signmat<<std::endl;
			temp=siteA.a[1]*signmat;
			kron(temp,siteB.c[1],temp2);
			temp2+=temp1;
			//std::cout<<temp2<<std::endl;
			//std::cout<<temp2.rowbase<<std::endl;
			//std::cout<<temp2.colbase<<std::endl;
			Transpose(temp2,temp);
			temp+=temp2;
			temp.scale(T);
		}
		Hinter.resize(temp.rowbase,temp.colbase);
		if(include_onsite=='y')
		{
			//U term
			temp1.resize(siteA.base,siteA.base);
			temp2.resize(siteB.base,siteB.base);
			if(sitenum<=2||addsiteposition=='l')
			{
				Mat_Mat_Mult(siteA.n[0],siteA.n[1],temp1);
				tempeye.geneye(siteB.base);
				kron(temp1,tempeye,Hinter);
			}
			if(sitenum<=2||addsiteposition=='r')
			{
				Mat_Mat_Mult(siteB.n[0],siteB.n[1],temp2);
				tempeye.geneye(siteA.base);
				kron(tempeye,temp2,temp1);
				Hinter+=temp1;
			}
			Hinter.scale(U);

			//std::cout<<Hinter<<std::endl;
			//Add
			Hinter+=temp;
			//for(int i=0;i<Hinter.rowbase.ordermap.size();i++)
			//    std::cout<<Hinter.rowbase.ordermap[i]<<std::endl;
		}
		else
		{
			Hinter=temp;
		}
	}
	else
	{
		std::cout<<"The evolution of hubbard chain is not supported yet."<<std::endl;
	}
}
 */

/*!
 \fn Block::calHinter_Heisenberg(Site &siteA,Site &siteB,char addsiteposition)
 */
/*
void Block::calHinter_Heisenberg(Site &siteA,Site &siteB,char addsiteposition)
{
	if(value_type=='r')
	{
		Rmatrix tempx,tempy;
		siteA.eval();
		siteB.eval();
		kron(siteA.c[0],siteB.a[0],tempx);
		//std::cout<<temp<<std::endl;
		kron(siteA.a[0],siteB.c[0],tempy);

		//std::cout<<temp2<<std::endl;
		tempx+=tempy;
		tempx.scale(T);
		//std::cout<<temp<<std::endl;
		if(addsiteposition=='r'&& sitenum==1)
		{
			std::cout<<"This is if code in Block::calHinter_Heisenberg is noly for temperary use to include the U term in impurity models, this must be replaced by more general method (matlab generated operators)"<<std::endl;
			Rmatrix temp1,temp2;
			temp1.geneye(siteA.base);
			temp2.geneye(siteB.base);
			temp1.scale(-0.5);
			temp2.scale(-0.5);
			temp1=siteA.n[0]+temp1;
			temp2=siteB.n[0]+temp2;
			kron(temp1,temp2,Hinter);
			Hinter.scale(1/3);
		}
		else
		{
			kron(siteA.n[0],siteB.n[0],Hinter);
			Hinter.scale(V);
		}
		Hinter+=tempx;
		// Hinter=tempx;
	}
	else
	{
		std::cout<<"calHinter_Heisenberg() has not been inplemented for complex numbers"<<std::endl;
	}
}
 */


