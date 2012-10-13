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


using namespace std;

double DTMat::MaxTruncateError=0;

DTMat::DTMat()
{
	value_type='r';
	KeptStatesNum=0;
}


DTMat::~DTMat()
{
}


DTMat::DTMat(char vtype,char hand)
{
	value_type=vtype;
	handside=hand;
}


/*!
\fn DTMat::DTMat(DTMat &mat)
 */
DTMat::DTMat(const DTMat &mat)
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


/*!
\fn DTMat::gendenmat(LaVectorDouble &wf)
 */
void DTMat::gendenmat(Rmatrix &wfmat,GQNBase &left,GQNBase &right)
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


void DTMat::gendenmat(Cmatrix &wfmat,GQNBase &left,GQNBase &right)
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


/*!
\fn DTMat::findtmat()
 */
void DTMat::findtmat(int tn)
{
	if(value_type=='r')
	{
		int Dim=denmatbase.Dim;
		KeptStatesNum=tn;
		if(Dim<=KeptStatesNum && MaxTruncateError==0)
		{
			//cout<<"!!!!!!!!!!!!!!!!!!!"<<endl;
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
		if(Dim<=KeptStatesNum && MaxTruncateError==0)
		{
			tmatbase=denmatbase;
			trunmatC.geneye(denmatbase);
			return;
		}
		int subnum=denmatbase.subnum;
		// cout<<"Make suer that Density Matrix is a block matrix"<<endl<<denmatC<<endl;
		LaGenMatComplex eigvec(Dim,Dim);
		LaVectorDouble eigval(Dim);

		eig(denmatC,eigvec,eigval);
        caltruncatemat(eigvec,eigval,trunmatC);

	}
}


/*!
    \fn DTMat::write(fout)
 */
void DTMat::write(ofstream &fout)
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


/*!
    \fn DTMat::read(fin)
 */
void DTMat::read(ifstream &fin)
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


/*!
    \fn DTMat::DTMat(ifstream &fin)
 */
DTMat::DTMat(ifstream &fin)
{
	read(fin);
}


/*!
    \fn DTMat::renorm()
 */
void DTMat::renorm()
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


/*!
    \fn DTMat::toComplex()
 */
void DTMat::toComplex()
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


/*!
    \fn DTMat::vonNeumannEntropy(LaVectorDouble &eigvals)
 */
double DTMat::vonNeumannEntropy(LaVectorDouble &eigvals)
{
	double entropy=0;
	int dim=eigvals.size();
	for(int i=0;i<dim;i++)
		if(eigvals(i)>1e-16)    entropy-=eigvals(i)*(log(eigvals(i))/log(2.0));

	double maxentropy=log(dim)/log(2.0);
	entropy=entropy/maxentropy;
	//if(entropy>0.5) cout<<eigvals<<endl;;
	return entropy;
}


/*!
    \fn DTMat::operator=(const DTMat& dtmat)
 */
DTMat& DTMat::operator=(const DTMat& dtmat)
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
