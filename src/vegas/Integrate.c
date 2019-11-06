/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Vegas
		last modified 13 Mar 15 th
*/


typedef struct {
  signature_t signature;
  count niter;
  number nsamples, neval;
  Cumulants cumul[];
} State;

static int Integrate3D(This *t, real *integral, real *error, real *prob)
{
  bin_t *bins;
  count dim, comp;
  int fail;

  StateDecl;
  csize_t statesize = sizeof(State) +
    NCOMP*sizeof(Cumulants) + NDIM*sizeof(Grid);
  Sized(State, state, statesize);
  Cumulants *c, *C = state->cumul + t->ncomp;
  Grid *state_grid = (Grid *)C;
  Array(Grid, margsum, NCOMP, NDIM);
  Vector(char, out, 128*NCOMP + 256);

  if( BadComponent(t) ) return -2;
  if( BadDimension(t) ) return -1;

  FrameAlloc(t, Master);
  ForkCores(t);
  Alloc(bins, t->nbatch*t->ndim);

  if( (fail = setjmp(t->abort)) ) goto abort;

  IniRandom(t);

  StateSetup(t);

  if( ini | ZAPSTATE ) {
    t->neval = 0;
    state->niter = 0;
    state->nsamples = t->nstart;
    FClear(state->cumul);
    if( ini ) GetGrid(t, state_grid);
  }

  /* main iteration loop */
  int cancel = 0;
  for( ; ; ) {
    number nsamples = state->nsamples;
    creal jacobian = 1./nsamples;

    FClear(margsum);

    for( ; nsamples > 0; nsamples -= t->nbatch ) {
      cnumber n = IMin(t->nbatch, nsamples);
      real *w = t->frame;
      real *x = w + n;
      real *f = x + n*t->ndim;
      real *lastf = f + n*t->ncomp;
      bin_t *bin = bins;

      while( x < f ) {
        real weight = jacobian;

        t->rng.getrandom(t, x);

        for( dim = 0; dim < t->ndim; ++dim ) {
          creal pos = *x*NBINS;
          ccount ipos = (count)pos;
          creal prev = (ipos == 0) ? 0 : state_grid[dim][ipos - 1];
          creal diff = state_grid[dim][ipos] - prev; 
          *x++ = prev + (pos - ipos)*diff;
          *bin++ = ipos;
          weight *= diff*NBINS;
        }

        *w++ = weight;
      }

     cancel = DoSample(t, n, w, f, t->frame, state->niter + 1);

     if (cancel)
     {
       fail = cancel;
       break;
     }

      bin = bins;
      w = t->frame;

      while( f < lastf ) {
        creal weight = *w++;
        Grid *m = &margsum[0][0];

        for( c = state->cumul; c < C; ++c ) {
          real wfun = weight*(*f++);
          if( wfun ) {
            c->sum += wfun;
            c->sqsum += wfun *= wfun;
            for( dim = 0; dim < t->ndim; ++dim )
              m[dim][bin[dim]] += wfun;
          }
          m += t->ndim;
        }

        bin += t->ndim;
      }
    }
    if (cancel) 
    { 
      fail = cancel;
      break; 
    }
    
    fail = 1;

    /* compute the integral and error values */

    for( c = state->cumul; c < C; ++c ) {
      real w = Weight(c->sum, c->sqsum, state->nsamples);
      real sigsq = 1/(c->weightsum += w);
      real avg = sigsq*(c->avgsum += w*c->sum);

      c->avg = LAST ? (sigsq = 1/w, c->sum) : avg;
      c->err = sqrtx(sigsq);
      //fail |= (c->err > MaxErr(c->avg));

      if( state->niter == 0 ) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w*c->sum;
      }
      c->chisq = c->chisqsum - avg*c->chisum;

      c->sum = c->sqsum = 0;
    }

    double* currentAverages = malloc((t->ncomp) * sizeof(double));

    for (int currentComp = 0; currentComp < t->ncomp; ++currentComp) {
      cCumulants *cc = &state->cumul[currentComp];
      currentAverages[currentComp] = cc->avg;
    }

    fail = t->updateFunc(currentAverages, t->ncomp, t->userdata);
    free(currentAverages);

    if (fail == 0 && t->neval >= t->mineval) break; 

    if( t->neval >= t->maxeval && !StateWriteTest(t) ) break;

    if( t->ncomp == 1 )
      for( dim = 0; dim < t->ndim; ++dim )
        RefineGrid(t, state_grid[dim], margsum[0][dim]);
    else {
      for( dim = 0; dim < t->ndim; ++dim ) {
        Grid wmargsum;
        Zap(wmargsum);
        for( comp = 0; comp < t->ncomp; ++comp ) {
          real w = state->cumul[comp].avg;
          if( w != 0 ) {
            creal *m = margsum[comp][dim];
            count bin;
            w = 1/Sq(w);
            for( bin = 0; bin < NBINS; ++bin )
              wmargsum[bin] += w*m[bin];
          }
        }
        RefineGrid(t, state_grid[dim], wmargsum);
      }
    }

    ++state->niter;
    state->nsamples += t->nincrease;

  }

  if (!cancel) 
  {
    for (comp = 0; comp < t->ncomp; ++comp) {
      cCumulants *c = &state->cumul[comp];
      integral[comp] = c->avg;
      error[comp] = c->err;
      prob[comp] = ChiSquare(c->chisq, state->niter);
    }
  }
  

abort:
  PutGrid(t, state_grid);
  free(bins);
  FrameFree(t, Master);

  StateRemove(t);

  return fail;
}

static int Integrate1D(This *t, real *integral, real *error, real *prob)
{
  bin_t *bins;
  count dim, comp;
  int fail;

  StateDecl;
  csize_t statesize = sizeof(State) +
    NCOMP * sizeof(Cumulants) + 1 * sizeof(Grid);
  Sized(State, state, statesize);
  Cumulants *c, *C = state->cumul + t->ncomp;
  Grid *state_grid = (Grid *)C;
  Array(Grid, margsum, NCOMP, 1);
  Vector(char, out, 128 * NCOMP + 256);

  if (BadComponent(t)) return -2;
  if (BadDimension(t)) return -1;

  FrameAlloc(t, Master);
  ForkCores(t);
  Alloc(bins, t->nbatch*t->ndim);

  if ((fail = setjmp(t->abort))) goto abort;

  IniRandom(t);

  StateSetup(t);

  if (ini | ZAPSTATE) {
    t->neval = 0;
    state->niter = 0;
    state->nsamples = t->nstart;
    FClear(state->cumul);
    if (ini) GetGrid(t, state_grid);
  }

  /* main iteration loop */
  int cancel = 0;
  for (; ; ) {
    number nsamples = state->nsamples;
    creal jacobian = 1. / nsamples;

    FClear(margsum);

    for (; nsamples > 0; nsamples -= t->nbatch) {
      cnumber n = IMin(t->nbatch, nsamples);
      real *w = t->frame;
      real *x = w + n;
      real *f = x + n * t->ndim;
      real *lastf = f + n * t->ncomp;
      bin_t *bin = bins;

      while (x < f) {
        real weight = jacobian;

        t->rng.getrandom(t, x);

        for (dim = 0; dim < t->ndim; ++dim) {
          creal pos = *x*NBINS;
          ccount ipos = (count)pos;
          creal prev = (ipos == 0) ? 0 : state_grid[dim][ipos - 1];
          creal diff = state_grid[dim][ipos] - prev;
          *x++ = prev + (pos - ipos)*diff;
          *bin++ = ipos;
          weight *= diff * NBINS;
        }

        *w++ = weight;
      }

      cancel = DoSample(t, n, w, f, t->frame, state->niter + 1);

      if (cancel)
      {
        fail = cancel;
        break;
      }

      bin = bins;
      w = t->frame;

      while (f < lastf) {
        creal weight = *w++;
        Grid *m = &margsum[0][0];

        for (c = state->cumul; c < C; ++c) {
          real wfun = weight * (*f++);
          if (wfun) {
            c->sum += wfun;
            c->sqsum += wfun *= wfun;
            for (dim = 0; dim < t->ndim; ++dim)
              m[dim][bin[dim]] += wfun;
          }
          m += t->ndim;
        }

        bin += t->ndim;
      }
    }
    if (cancel)
    {
      fail = cancel;
      break;
    }

    fail = 1;

    /* compute the integral and error values */

    for (c = state->cumul; c < C; ++c) {
      real w = Weight(c->sum, c->sqsum, state->nsamples);
      real sigsq = 1 / (c->weightsum += w);
      real avg = sigsq * (c->avgsum += w * c->sum);

      c->avg = LAST ? (sigsq = 1 / w, c->sum) : avg;
      c->err = sqrtx(sigsq);
      //fail |= (c->err > MaxErr(c->avg));

      if (state->niter == 0) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w * c->sum;
      }
      c->chisq = c->chisqsum - avg * c->chisum;

      c->sum = c->sqsum = 0;
    }

    double* currentAverages = malloc((t->ncomp) * sizeof(double));

    for (int currentComp = 0; currentComp < t->ncomp; ++currentComp) {
      cCumulants *cc = &state->cumul[currentComp];
      currentAverages[currentComp] = cc->avg;
    }

    fail = t->updateFunc(currentAverages, t->ncomp, t->userdata);
    free(currentAverages);

    if (fail == 0 && t->neval >= t->mineval) break;

    if (t->neval >= t->maxeval && !StateWriteTest(t)) break;

    if (t->ncomp == 1)
      for (dim = 0; dim < t->ndim; ++dim)
        RefineGrid(t, state_grid[dim], margsum[0][dim]);
    else {
      for (dim = 0; dim < t->ndim; ++dim) {
        Grid wmargsum;
        Zap(wmargsum);
        for (comp = 0; comp < t->ncomp; ++comp) {
          real w = state->cumul[comp].avg;
          if (w != 0) {
            creal *m = margsum[comp][dim];
            count bin;
            w = 1 / Sq(w);
            for (bin = 0; bin < NBINS; ++bin)
              wmargsum[bin] += w * m[bin];
          }
        }
        RefineGrid(t, state_grid[dim], wmargsum);
      }
    }

    ++state->niter;
    state->nsamples += t->nincrease;

  }

  if (!cancel)
  {
    for (comp = 0; comp < t->ncomp; ++comp) {
      cCumulants *c = &state->cumul[comp];
      integral[comp] = c->avg;
      error[comp] = c->err;
      prob[comp] = ChiSquare(c->chisq, state->niter);
    }
  }


abort:
  PutGrid(t, state_grid);
  free(bins);
  FrameFree(t, Master);

  StateRemove(t);

  return fail;
}

