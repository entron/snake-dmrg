#include "dtmat.h"
#include "public.h"

#include <laindex.h>
#include <blas1pp.h>
#include<blas3pp.h>
#include<blaspp.h>
#include<laslv.h>

#include <iostream>
#include<cstring>
#include<algorithm>


snake::physics::DTMat::DTMat()
{
	value_type='r';
	KeptStatesNum=0;
}


snake::physics::DTMat::~DTMat()
{
}


snake::physics::DTMat::DTMat(char vtype,char hand)
{
	value_type=vtype;
	handside=hand;
}



snake::physics::DTMat::DTMat(const DTMat &mat)
{
	denmat=mat.denmat;
	denmatC=mat.denmatC;
	denmatbase=mat.denmatbase;
	handside=mat.handside;
	leftbase=mat.leftbase;
	rightbase=mat.rightbase;
	tmatbase=mat.tmatbase;
	trunmat=mat.trunmat;
	trunmatC=mat.trunmatC;
	//trunmattrans=mat.trunmattrans;
	//trunmattrans2=mat.trunmattrans2;
	KeptStatesNum=mat.KeptStatesNum;
	value_type=mat.value_type;
}



void snake::physics::DTMat::gendenmat(Rmatrix &wfmat, snake::math::GQNBase &left, snake::math::GQNBase &right)
{
	leftbase=left;
	rightbase=right;

	if(handside=='l')
	{
		denmatbase=leftbase;
		denmat.resize(denmatbase,denmatbase);
		Mat_Trans_Mat_Mult(wfmat,wfmat,denmat);
	}
	else
	{
		denmatbase=rightbase;
		denmat.resize(denmatbase,denmatbase);
		Mat_Mat_Trans_Mult(wfmat,wfmat,denmat);
	}

}


void snake::physics::DTMat::gendenmat(Cmatrix &wfmat, snake::math::GQNBase &left, snake::math::GQNBase &right)
{
	leftbase=left;
	rightbase=right;

	if(handside=='l')
	{
		denmatbase=leftbase;
		denmatC.resize(denmatbase,denmatbase);
		Mat_Trans_Mat_Mult(wfmat,wfmat,denmatC);
	}
	else
	{
		denmatbase=rightbase;
		denmatC.resize(denmatbase,denmatbase);
		Mat_Mat_Trans_Mult(wfmat,wfmat,denmatC);
	}
}


void snake::physics::DTMat::findtmat(int tn)
{
	if(value_type=='r')
	{
		int Dim=denmatbase.Dim;
		KeptStatesNum=tn;
        if(Dim<=KeptStatesNum && Max_Truncate_Error==0)
		{
			//std::cout<<"!!!!!!!!!!!!!!!!!!!"<<std::endl;
			tmatbase=denmatbase;
			trunmat.geneye(denmatbase);
			return;
		}
		int subnum=denmatbase.subnum;
		LaGenMatDouble eigvec(Dim,Dim);
		LaVectorDouble eigval(Dim);

		eig(denmat,eigvec,eigval);


        caltruncatemat(eigvec,eigval,trunmat);

	}
	else
	{
		int Dim=denmatbase.Dim;
		KeptStatesNum=tn;
        if(Dim<=KeptStatesNum && Max_Truncate_Error==0)
		{
			tmatbase=denmatbase;
			trunmatC.geneye(denmatbase);
			return;
		}
		int subnum=denmatbase.subnum;
		// std::cout<<"Make suer that Density Matrix is a block matrix"<<std::endl<<denmatC<<std::endl;
		LaGenMatComplex eigvec(Dim,Dim);
		LaVectorDouble eigval(Dim);

		eig(denmatC,eigvec,eigval);
        caltruncatemat(eigvec,eigval,trunmatC);

	}
}



void snake::physics::DTMat::write(std::ofstream &fout)
{

	fout.write(&value_type,sizeof value_type);
	fout.write(&handside,sizeof handside);
	leftbase.write(fout);
	rightbase.write(fout);
	tmatbase.write(fout);
	if(value_type=='r')
	{
		trunmat.write(fout);
		//trunmattrans.write(fout);
	}
	else
	{
		trunmatC.write(fout);
		//trunmattrans2.write(fout);
	}
}



void snake::physics::DTMat::read(std::ifstream &fin)
{
	fin.read(&value_type,sizeof value_type);
	fin.read(&handside,sizeof handside);
	leftbase.read(fin);
	rightbase.read(fin);
	tmatbase.read(fin);
	if(value_type=='r')
	{
		trunmat.read(fin);
		//trunmattrans.read(fin);
	}
	else
	{
		trunmatC.read(fin);
		//trunmattrans2.read(fin);
	}

}



snake::physics::DTMat::DTMat(std::ifstream &fin)
{
	read(fin);
}



void snake::physics::DTMat::renorm()
{
	if(value_type=='r')
	{
		Rmatrix temp(tmatbase,denmatbase);
		Mat_Trans_Mat_Mult(trunmat,denmat,temp);
		denmat.resize(tmatbase,tmatbase);
		Mat_Mat_Mult(temp,trunmat,denmat);
		denmatbase=tmatbase;
	}
	else
	{
		Cmatrix temp(tmatbase,denmatbase);
		Mat_Trans_Mat_Mult(trunmatC,denmatC,temp);
		denmatC.resize(tmatbase,tmatbase);
		Mat_Mat_Mult(temp,trunmatC,denmatC);
		denmatbase=tmatbase;
	}

}



void snake::physics::DTMat::toComplex()
{
	if(value_type=='r')
	{
		value_type='c';
		R2C(denmat,denmatC);
		R2C(trunmat,trunmatC);
		//R2C(trunmattrans,trunmattrans2);
		denmat.delmat();
		trunmat.delmat();
		//trunmattrans.delmat();
	}
}



double snake::physics::DTMat::vonNeumannEntropy(LaVectorDouble &eigvals)
{
	double entropy=0;
	int dim=eigvals.size();
	for(int i=0;i<dim;i++)
		if(eigvals(i)>1e-16)    entropy-=eigvals(i)*(log(eigvals(i))/log(2.0));

	double maxentropy=log(dim)/log(2.0);
	entropy=entropy/maxentropy;
	//if(entropy>0.5) std::cout<<eigvals<<std::endl;;
	return entropy;
}



snake::physics::DTMat& snake::physics::DTMat::operator=(const DTMat& dtmat)
{
	value_type=dtmat.value_type;
	denmat=dtmat.denmat;
	denmatC=dtmat.denmatC;
	denmatbase=dtmat.denmatbase;
	handside=dtmat.handside;
	leftbase=dtmat.leftbase;
	rightbase=dtmat.rightbase;
	tmatbase=dtmat.tmatbase;
	trunmat=dtmat.trunmat;
	//trunmattrans=dtmat.trunmattrans;
	trunmatC=dtmat.trunmatC;
	//trunmattrans2=dtmat.trunmattrans2;
	KeptStatesNum=dtmat.KeptStatesNum;
}
