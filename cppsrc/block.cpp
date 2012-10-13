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
Block::Block(string &filename)
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
 cout<<"!!!!!!!!!!!!!!!!Be Careful about the potential on the impurity!!!!!!!!!!!!!!!!!!!!!!"<<endl;
 cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Edit block.cpp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
 cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!This is a alpha version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
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
 cout<<"Block dim is "<<base.Dim<<endl;
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
 cout<<"Block dim is "<<base.Dim<<endl;
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
 cout<<"Block dim is "<<base.Dim<<endl;
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
 //cout<<"Hinter"<<Hinter<<endl;

 }
 */
void Block::calHinter(Site &siteA,Site &siteB,char addsiteposition,double HoppingT, double OnSiteE, double TwoSitesV)
{
	//cout<<"HoppingT="<<HoppingT<<endl;


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
			cout<<"OnsiteE="<<OnSiteE<<endl;
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
  	//cout<<Hinter<<endl;

		//Two sites interaction term
		//cout<<"!!!!!!!!!!Model Dependent Part!!!!!!!!!!"<<endl;
		//cout<<"H_i=V(n_i-1/2)(n_{i+1}-1/2)"<<endl;
		if(TwoSitesV!=0)
		{
			cout<<"TwoSitesV="<<TwoSitesV<<endl;
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
		cout<<"calHinter_Heisenberg() has not been inplemented for complex numbers"<<endl;
	}
  //cout<<Hinter<<endl;
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
	string filename;
	stringstream stemp;
	stemp<<sitenum;
	filename=prefix+stemp.str();
  //cout<<endl<<filename<<endl;
	ofstream fout(filename.c_str(),ios_base::out|ios_base::binary);
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
void Block::read(string &filename)
{
	ifstream fin(filename.c_str(),ios_base::in|ios_base::binary);
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
 \fn Block::calN(ofstream fout)
 */
void Block::calN(ofstream &fout,char hand)
{
	if(value_type=='r')
	{
		Rmatrix tempmat(base,base);
		dtmat->renorm();
		if(hand=='l')
			for(int i=0;i<sitenum;i++)
			{
				//cout<<"Calculate Site "<<i<<endl;
				//cout<<dtmat->denmat.size(0)<<endl;
				//cout<<site[i].n.size(0)<<endl;
				for(int j=0;j<site[i].num;j++)
				{
					Mat_Mat_Mult(dtmat->denmat,site[i].n[j],tempmat);
					fout<<tempmat.trace()<<" ";
				}
				fout<<endl;
				//cout<<trace(tempmat)<<endl;
			}
		else
			for(int i=sitenum-1;i>=0;i--)
			{
				for(int j=0;j<site[i].num;j++)
				{
					Mat_Mat_Mult(dtmat->denmat,site[i].n[j],tempmat);
					fout<<tempmat.trace()<<" ";
				}
				fout<<endl;
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
				//cout<<"Calculate Site "<<i<<endl;
				//cout<<dtmat->denmat.size(0)<<endl;
				//cout<<site[i].n.size(0)<<endl;
				for(int j=0;j<site[i].num;j++)
				{
					Mat_Mat_Mult(dtmat->denmatC,site[i].nC[j],tempmat2);
					C2R(tempmat2,tempmat);
					fout<<tempmat.trace()<<" ";
				}
				cout<<endl;
				//cout<<trace(tempmat)<<endl;
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
				cout<<endl;
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
			//cout<<siteA.a[0]<<endl;
			//cout<<siteB.c[0]<<endl;
			kron(siteA.a[0],siteB.c[0],temp1);
			//cout<<temp1<<endl;
			kron(siteA.a[1],siteB.c[1],temp);
			//cout<<temp<<endl;
			temp+=temp1;
			//cout<<temp<<endl;
			//cout<<temp.rowbase<<endl;
			//cout<<temp.colbase<<endl;
			Transpose(temp,temp2);
			//cout<<temp<<endl;
			//cout<<temp2<<endl;
			temp+=temp2;
			temp.scale(T);
			//cout<<temp<<endl;
		}
		else
		{
			//Hopping term
			Rmatrix signmat;
			signmat.gensignmat(siteA.base,0);
			//cout<<signmat<<endl;
			temp=siteA.a[0]*signmat;
			kron(temp,siteB.c[0],temp1);
			signmat.gensignmat(siteA.base,1);
			//cout<<signmat<<endl;
			temp=siteA.a[1]*signmat;
			kron(temp,siteB.c[1],temp2);
			temp2+=temp1;
			//cout<<temp2<<endl;
			//cout<<temp2.rowbase<<endl;
			//cout<<temp2.colbase<<endl;
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

			//cout<<Hinter<<endl;
			//Add
			Hinter+=temp;
			//for(int i=0;i<Hinter.rowbase.ordermap.size();i++)
			//    cout<<Hinter.rowbase.ordermap[i]<<endl;
		}
		else
		{
			Hinter=temp;
		}
	}
	else
	{
		cout<<"The evolution of hubbard chain is not supported yet."<<endl;
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
		//cout<<temp<<endl;
		kron(siteA.a[0],siteB.c[0],tempy);

		//cout<<temp2<<endl;
		tempx+=tempy;
		tempx.scale(T);
		//cout<<temp<<endl;
		if(addsiteposition=='r'&& sitenum==1)
		{
			cout<<"This is if code in Block::calHinter_Heisenberg is noly for temperary use to include the U term in impurity models, this must be replaced by more general method (matlab generated operators)"<<endl;
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
		cout<<"calHinter_Heisenberg() has not been inplemented for complex numbers"<<endl;
	}
}
 */


