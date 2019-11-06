/*
	Vegas.c
		Vegas Monte Carlo integration
		by Thomas Hahn
		last modified 25 Nov 14 th
*/


#define VEGAS
#define ROUTINE "Vegas"

#include "decl.h"
#include "CSample.c"

#include <cuba_export.h>

/*********************************************************************/

Extern CUBA_EXPORT void(Vegas)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata, CustomSampleFunc customSample,UpdateFunc uf, cnumber nvec,
  creal epsrel, creal epsabs, cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  cnumber nbatch, cint gridno,
  cchar *statefile, Spin **pspin,
  number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  This t;
  
  t.sampleFunc = customSample;
  t.updateFunc = uf;
  t.ndim = ndim;
  t.ncomp = ncomp;
  t.integrand = integrand;
  t.userdata = userdata;
  t.nvec = nvec;
  t.epsrel = epsrel;
  t.epsabs = epsabs;
  t.flags = flags;
  t.seed = seed;
  t.mineval = mineval;
  t.maxeval = maxeval;
  t.nstart = nstart;
  t.nincrease = nincrease;
  t.nbatch = nbatch;
  t.gridno = gridno;
  t.statefile = statefile;
  if (ndim == 1) {
    *pfail = Integrate1D(&t, integral, error, prob);
  }
  else {
    *pfail = Integrate3D(&t, integral, error, prob);
  }
  *pneval = t.neval;

}

