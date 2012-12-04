#ifndef SUPBLOCK_H
#define SUPBLOCK_H

/**
 This class contains information related to the super block.
 This class might be rewrited so that it's inhereted from Block class

 @author Cheng Guo
 */


#include<fstream>
#include<lavd.h>
#include<vector>
#define LA_COMPLEX_SUPPORT
#include<lavc.h>
#include <complex>


#include "gqnbase.h"
#include "setting.h"
#include "gqnmat.h"
#include "gqn.h"
#include"block.h"
#include "dtmat.h"
#include "blocham.h"
#include "site.h"

using namespace snake::math;

namespace snake
{

namespace physics
{
extern long int multnum;
extern std::vector<snake::physics::Site> allfreesites;

class SupBlock {
public:
    char value_type;
    int KeptStatesNum;
    std::vector<snake::math::GQN> TargetGQN;
    std::vector<snake::math::GQN> TargetGQN2;
    int sitenum;

    ///Left and right block.
    Block *L,*R,*oldL,*oldR;

    ///Ground state wave function.
    LaVectorDouble wf;
    LaVectorComplex wfC;
    LaVectorComplex wfC2;
    //GQN tnum2;
    ///Left and right blocks' dtmats,useful when time evolve
    std::vector<DTMat> leftdtmat,rightdtmat;

    ///Matrix Form of wavefunction
    Rmatrix wfmat;
    Cmatrix wfmatC;
    Cmatrix wfmatC2;

    snake::math::GQNBase oldleftbase,oldrightbase,leftbase,rightbase;

private:


    int *mapdim;
    int **map;
    int ****index;
    int H2Dim;
    Site midsite1, midsite2,freesite;

    int mtimes;
    LaGenMatDouble Hlr;
    std::ofstream fout_1stsite_n_t, fout_entropy_t, fout_steperror_t;
public:
    SupBlock();
    ~SupBlock();
    //SupBlock(Block *left,Block *right,Block *oleft,Block *oright);
    ///For interaction Terms are different
    SupBlock(Block *left,Block *right,Block *oleft,Block *oright,LaGenMatDouble &Hi);
    ///Find the ground state of target good quantum number tTargetGQN.
    void CalGroundState();

    ///Renormalize left and right block.
    void renorm(int KeptStatesNum);

    ///Renormalize right block.
    void renormright(int KeptStatesNum);
    void renormleft(int KeptStatesNum);


    // void calN(char *filename);

    ///Calculate correlation fucntion
    void calCF(char* filename);

    ///Wave function transformation using White's method
    void moveleft(DTMat &leftdtmat,DTMat &rightdtmat);
    void moveright(DTMat &leftdtmat,DTMat &rightdtmat);

    ///Time (real or imaginary) evolve
    void evolve(std::vector<LaGenMatDouble> &UT,int timesteps);
    void evolve(std::vector<LaGenMatComplex> &Ut, int timesteps);
    void evolve(std::vector<LaGenMatComplex> &ut, std::vector<LaGenMatComplex> &Ut, int timesteps);
    /**Read data from files wrote by FINITE DMRG etc. so that supblock
     *can time evolve. n is the number of site which are exact.*/
    void prepare();

    ///Read dtmats from file saved during the finite DMRG step
    void loaddtmats();

    void toComplex();

    /**Sweep from middle to left and do nothing, useful when evolve
     *the finite DMRG wf */
    void sweep2left(int);

    void sweep2leftmost();

    ///The middle two free site operator T multiply w
    ///out=out+TO*in
    template<class MATType, class Type>
    void middlemult(MATType &TO,Type *in,Type *out)
    {
        //std::cout<<TO<<std::endl;
        for (int i=0;i<H2Dim;i++)
            for (int k=0;k<H2Dim;k++)
                if (TO(i,k) != 0)
                    for (int j=0;j<mapdim[i];j++)
                        //out[map[i][j]]+=TO(i,k)*in[map[k][j]];
                        out[map[i][j]]=out[map[i][j]]+TO(i,k)*in[map[k][j]]; //For the genrality to sacrify a very little bit of efficiency.
    }

    template<class VType,class MATType>
    void middlemult(MATType &TO,VType &f)
    {
        VType out;
        out=f;
        f=0;
        middlemult(TO,out.addr(),f.addr());
        multnum++;
    }


    ///Apply local site operator on the wf
    void applyop(LaGenMatComplex &op,int thesite);
    void write(char *filename);
    void read(char *filename);

    ///Generate mapdim and map
    void genmiddlemap(std::vector<snake::math::GQN> &tgqn);

    ///Generate two freesite index.
    void genindex();
    void delindex();
    //  void applyOPonDot(Rmatrix &OP);
    void creatoutputfiles();
    void closeoutputfiles();

    /*!
     \fn SupBlock::evalwavefuncmat()
     */
    template<class WFType,class WFMATType>
    void evalwfmat(WFType &f,WFMATType &mat, std::vector<snake::math::GQN> TGQN)
    {
        int vstart=0;
        //std::cout<<f.size()<<std::endl;
        //if(mat.rowbase!=rightbase || mat.colbase!=leftbase)
        mat.resize(rightbase,leftbase);
        //std::cout<<rightbase<<std::endl;
        //std::cout<<leftbase<<std::endl;
        //std::cout<<f<<std::endl;
        WFType tempv;
        typename WFMATType::value_type tempmat;
        for (int i=0;i<leftbase.subnum;i++)
        {
            for (int j=0;j<rightbase.subnum;j++)
            {
                if (leftbase.subgqn[i]+rightbase.subgqn[j]==TGQN)
                {
                    tempv=f(LaIndex(vstart,vstart+leftbase.dim[i]*rightbase.dim[j]-1));
                    vstart+=leftbase.dim[i]*rightbase.dim[j];
                    tempmat.resize(rightbase.dim[j],leftbase.dim[i]);
                    snake::math::vecCmat(tempv,tempmat);
                    //std::cout<<tempmat<<std::endl;
                    mat.subnum++;
                    mat.submat.resize(mat.subnum);
                    mat.pmat(j,i)=mat.subnum-1;
                    mat.submat[mat.subnum-1]=tempmat;
                }
            }
        }
    }


    /*!
     \fn SupBlock::extractwf(int TargetGQN)
     */
    template<class WFType,class WFMATType>
    void extractwf(WFMATType &mat,WFType &f, std::vector<snake::math::GQN> TGQN)
    {
        WFType tempv;
        typename WFMATType::value_type tempmat;
        f.resize(0,1);
        //printvector(TargetGQN);std::cout<<std::endl;
        for (int i=0;i<leftbase.subnum;i++)
        {
            for (int j=0;j<rightbase.subnum;j++)
            {
                if (leftbase.subgqn[i]+rightbase.subgqn[j]==TGQN)
                {
                    if (mat.pmat(j,i)!=-1)
                        snake::math::mat2vec(tempv,mat.submat[mat.pmat(j,i)]);
                    else
                    {
                        tempv.resize(leftbase.dim[i]*rightbase.dim[j],1);
                        tempv=0;
                    }
                    snake::math::join(f,tempv);
                }
            }
        }
        //std::cout<<"Wave func is "<<std::endl<<wf<<std::endl;
    }

private:
    ///This function is used move a middle site form one side block to the other
    template<class Type>
    void reshapewfmat(Type &wfmatfull,char hand)
    {
        //std::cout<<wfmat<<std::endl;
        int matrow=wfmatfull.size(0);
        int matcol=wfmatfull.size(1);
        int sitedim=freesite.base.Dim;
        int newrow, newcol;
        Type temp,tempsub;
        if (hand=='r')
        {
            newrow=matrow/sitedim;
            newcol=matcol*sitedim;
            temp.resize(newrow,newcol);
            int col=0;
            for (int i=0;i<matcol;i++)
                for (int j=0;j<sitedim;j++)
                {
                    tempsub=wfmatfull(LaIndex(j*newrow,j*newrow+newrow-1),LaIndex(i));
                    //std::cout<<tempsub<<std::endl;
                    temp(LaIndex(),LaIndex(col)).inject(tempsub);
                    col++;
                }
        }
        else
        {
            newrow=matrow*sitedim;
            newcol=matcol/sitedim;
            temp.resize(newrow,newcol);
            int col=0;
            for (int i=0;i<newcol;i++)
                for (int j=0;j<sitedim;j++)
                {
                    tempsub=wfmatfull(LaIndex(),LaIndex(col));
                    temp(LaIndex(j*matrow,j*matrow+matrow-1),LaIndex(i)).inject(tempsub);
                    col++;
                }
        }
        wfmatfull=temp;
        // std::cout<<wfmat<<std::endl;
    }



    ///Renorm one base of wavefunction matrix
    void renormwfmat(DTMat &dtmat);
    void unrenormwfmat(DTMat &dtmat);



    /**The following three functions is to setup the temperature evolve
     *initial conditions*/
    void geninitialwf();
    void geninitialbases();
    void geninitialdtmats(std::vector<DTMat> &left, std::vector<DTMat> &right);

    ///Correlation function.<f|m1 m0|f>. hand tells which side m belongs.
    double corrfunc(Rmatrix &m1,char hand1,Rmatrix &m0,char hand0);

    void av(int n,double *in,double *out);

    /**This function utilize arpack to calculate the eigensystem of super
     *block hamiltonian.
     *I copy dsaupd funciton from:
     *http://cyclone.harvard.edu/people/shaw/programs/lapack.html
     */
    void dsaupd(int n,int nev,double *Evals,double **Evecs);
    int calDim();
    void leftmult(double *in ,double *out);
    void rightmult(double *in,double *out);
    void deletemiddlemap();
    void onetimestep(LaGenMatComplex& impurity_OP_t, std::vector<LaGenMatComplex>& rt_OP);
};

extern "C" void dsaupd_(int *ido, char *bmat, int *n, char *which,
                            int *nev, double *tol, double *resid, int *ncv,
                            double *v, int *ldv, int *iparam, int *ipntr,
                            double *workd, double *workl, int *lworkl,
                            int *info);

extern "C" void dseupd_(int *rvec, char *All, int *select, double *d,
                            double *v, int *ldv0, double *sigma,
                            char *bmat, int *n, char *which, int *nev,
                            double *tol, double *resid, int *ncv, double *v1,
                            int *ldv, int *iparam, int *ipntr, double *workd,
                            double *workl, int *lworkl, int *ierr);

extern "C"
{
    void dnaupd_(int *ido, const char *bmat, const int *n, char *which,const int *nev, double *tol, double *resid,const int *ncv, double *v, const int *ldv,int *iparam, int *ipntr, double *workd,double *workl, const int *lworkl, int *info);

    void dneupd_(int *rvec, char *howmny, int *select, double *dr, double *di, double *v0, int *ldv0, double *sigmar,double *sigmai, double *workev, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);
}


//extern std::ofstream foccnum,fentropy;
}}

#endif
