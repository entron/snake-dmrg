#ifndef DTMAT_H
#define DTMAT_H

/**
Density and truncation matrix class.The main purpose is to evaluate density matrix,and truncation matrix.
@author Cheng Guo
 */

#include "gqnmat.h"
#include "setting.h"
#include "gqnbase.h"

#include<gmd.h>
#include <lavd.h>
#define LA_COMPLEX_SUPPORT
#include <gmc.h>

#include<vector>
#include<fstream>

namespace snake
{

namespace physics
{

class DTMat{
public:
	Rmatrix denmat,trunmat;
	Cmatrix denmatC,trunmatC;

	///Truncation number.
	int KeptStatesNum;

	char handside;

	char value_type;

    snake::math::GQNBase leftbase,rightbase,denmatbase,tmatbase;
	double entropy;

public:
	DTMat();
	~DTMat();
	DTMat(char vtype,char hand);
    DTMat(std::ifstream &fin);
	DTMat(const DTMat &mat);

    void gendenmat(Rmatrix &wfmat,snake::math::GQNBase &l,snake::math::GQNBase &r);
    void gendenmat(Cmatrix &wfmat,snake::math::GQNBase &l,snake::math::GQNBase &r);

	void findtmat(int tn);

    void write(std::ofstream &fout);
    void read(std::ifstream &fin);

	void renorm();

	void toComplex();
	///Note that eigvals should be normalized first
	double vonNeumannEntropy(LaVectorDouble &eigvals);
	DTMat& operator=(const DTMat& dtmat);
private:
	template<typename DMType, typename MType>
	void eig(DMType &dm, MType &evc, LaVectorDouble &evl);

	/**Find cutedge,which is the quantity to justify whether to keep an
	 *eigenvector of density matrix.
	 */
	template<typename TMType, typename MType>
	void caltruncatemat(MType &vec, LaVectorDouble &evl, TMType &trumat);

};


template<typename DMType, typename MType>
void snake::physics::DTMat::eig(DMType &dm, MType &evc, LaVectorDouble &evl)
{
	//evc.resize(Dim,Dim);
	//evl.resize(Dim);
	evl=0;
	evc=0;
	int subnum=denmatbase.subnum;

	int start=0;
	for(int i=0;i<subnum;i++)
	{
		int subdim=denmatbase.dim[i];
		MType subdenmat(subdim,subdim),subeigvec(subdim,subdim);
		LaVectorDouble subeigval(subdim);
		if(dm.pmat(i,i)==-1)
		{
			subdenmat.resize(denmatbase.dim[i],denmatbase.dim[i]);
			subdenmat=0;
		}
		else
			subdenmat.copy(dm.submat[dm.pmat(i,i)]);
        snake::math::SSMED(subdenmat.addr(),subdim,subeigval.addr());
		evl(LaIndex(start,start+subdim-1)).inject(subeigval);
		evc(LaIndex(start,start+subdim-1),LaIndex(start,start+subdim-1)).inject(subdenmat);
		start+=subdim;
	}

	double mod=Blas_Norm1(evl);
	evl.scale(1/mod);
}

template<typename TMType, typename MType>
void snake::physics::DTMat::caltruncatemat(MType &evc, LaVectorDouble &evl, TMType &trumat)
{
	int Dim=denmatbase.Dim;

  //std::cout<<evl<<std::endl;
    std::vector<double> temp(Dim);
	for(int i=0;i<Dim;i++)
		temp[i]=evl(i);
	sort(temp.begin(),temp.end());

	double cutedge;

    if(Dim<=KeptStatesNum && Max_Truncate_Error>0)
        cutedge=Max_Truncate_Error;
	else
	{
		cutedge=temp[Dim-KeptStatesNum];
        if(cutedge<Max_Truncate_Error)
            cutedge=Max_Truncate_Error;
	}

	double trunerror=0;
	for(int i=0;temp[i]<cutedge;i++)
		trunerror+=temp[i];
	entropy=vonNeumannEntropy(evl);

	///Generate trunmat
	int *mark;
	mark=new int [denmatbase.Dim];
	for(int i=0;i<denmatbase.Dim;i++) mark[i]=0;
	truncate(denmatbase,tmatbase,evl,cutedge,mark);
    //for(int i=0;i<denmatbase.Dim;i++) std::cout<<mark[i]<<" ";
    //std::cout<<std::endl

	MType trunmatfull(denmatbase.Dim,tmatbase.Dim);
	int q=0;
	for(int i=0;i<denmatbase.Dim;i++)
		if(mark[i]==1)
		{
			trunmatfull.col(q).inject(evc.col(i));
			q++;
		}
	delete [] mark;

    //std::cout<<trunmatfull<<std::endl;
    //std::cout<<denmatbase<<std::endl;
    //std::cout<<tmatbase<<std::endl;
	trumat=TMType(trunmatfull,denmatbase,tmatbase,1);

}

}}



#endif
