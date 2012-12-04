#include "gqnbase.h"
#include "setting.h"


void snake::math::R2C(snake::math::GQNMat<LaGenMatDouble> &r,snake::math::GQNMat<LaGenMatComplex> &c)
{
  c.resize(r.rowbase,r.colbase);
  c.pmat=r.pmat;
  c.subnum=r.subnum;
  c.submat.resize(c.subnum);
  for(int i=0;i<c.subnum;i++)
    c.submat[i]=r.submat[i].to_LaGenMatComplex();
}

void snake::math::C2R(snake::math::GQNMat<LaGenMatComplex> &c,snake::math::GQNMat<LaGenMatDouble> &r)
{
  r.resize(c.rowbase,c.colbase);
  r.pmat=c.pmat;
  r.subnum=c.subnum;
  r.submat.resize(r.subnum);
  for(int i=0;i<r.subnum;i++)
    r.submat[i]=c.submat[i].real_to_LaGenMatDouble();
}
