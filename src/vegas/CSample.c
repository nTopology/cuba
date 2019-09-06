/*
	CSample.c
		the serial sampling routine
		for the C versions of the Cuba routines
		by Thomas Hahn
		last modified 9 Oct 14 th
*/


static coreinit cubafun_;
extern int cubaverb_;
extern corespec cubaworkers_;


static inline number SampleRaw(This *t, number n, creal *x, real *f,
  cint core VES_ONLY(, creal *w, ccount iter))
{
  return t->sampleFunc(n, x, f, core, w, iter, t->integrand, &t->ndim, &t->ncomp, t->userdata);
 
}

/*********************************************************************/

static inline int DoSampleSerial(This *t, cnumber n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter))
{
  MasterInit();
  t->neval += n;
  return SampleRaw(t, n, x, f, 0 VES_ONLY(, w, iter));
}

/*********************************************************************/

#ifdef HAVE_FORK

static void DoSample(This *t, number n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter));
DIV_ONLY(static int Explore(This *t, cint iregion);)

#else

#define DoSample DoSampleSerial
#define Explore ExploreSerial
#define ForkCores(t)

static inline void WaitCores(This *t, Spin **pspin)
{
  if( Invalid(pspin) ) MasterExit();
}

#define WaitCores(t, pspin)

#endif

#ifdef DIVONNE
static inline count SampleExtra(This *t, cBounds *b)
{
  number n = t->nextra;
  t->peakfinder(&t->ndim, b, &n, t->xextra, t->userdata);
  DoSample(t, n, t->xextra, t->fextra);
  return n;
}
#endif

#include "common.c"

#ifdef HAVE_FORK
#include "Parallel.c"
#endif

#include "Integrate.c"

