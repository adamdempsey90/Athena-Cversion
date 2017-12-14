#include "copyright.h"
/*============================================================================*/
/*! \file convection.c
 *  \brief Problem generator for RT instabilty.
 *
 * PURPOSE: Problem generator for RT instabilty.  Gravitational pot. is
 *   hardwired to be 0.1z. Density difference is hardwired to be 2.0 in 2D, and
 *   is set by the input parameter <problem>/rhoh in 3D (default value is 3.0).
 *   This reproduces 2D results of Liska & Wendroff, 3D results of
 *   Dimonte et al.
 * 
 * FOR 2D HYDRO:
 * Problem domain should be -1/6 < x < 1/6; -0.5 < y < 0.5 with gamma=1.4 to
 * match Liska & Wendroff. Interface is at y=0; perturbation added to Vy
 * Gravity acts in the y-direction.  Special reflecting boundary conditions
 *   added in x2 to improve hydrostatic eqm (prevents launching of weak waves)
 * Atwood number A = (d2-d1)/(d2+d1) = 1/3
 *
 * FOR 3D:
 * Problem domain should be -.05 < x < .05; -.05 < y < .05, -.1 < z < .1
 * Use gamma=5/3 to match Dimonte et al.
 * Interface is at z=0; perturbation added to Vz
 * Gravity acts in the z-direction.  Special reflecting boundary conditions
 *   added in x3 to improve hydrostatic eqm (prevents launching of weak waves)
 * Atwood number A = (d2-d1)/(d2+d1) = 1/2
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - ran2() - random number generator from NR
 * - reflect_ix2() - sets BCs on L-x2 (left edge) of grid used in 2D
 * - reflect_ox2() - sets BCs on R-x2 (right edge) of grid used in 2D
 * - reflect_ix3() - sets BCs on L-x3 (left edge) of grid used in 3D
 * - reflect_ox3() - sets BCs on R-x3 (right edge) of grid used in 3D
 * - grav_pot2() - gravitational potential for 2D problem (accn in Y)
 * - grav_pot3() - gravitational potential for 3D problem (accn in Z)
 *
 * REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)    */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2() - random number generator from NR
 * flux_ix2() - sets BCs on L-x2 (left edge) of grid used in 2D
 * flux_ox2() - sets BCs on R-x2 (right edge) of grid used in 2D
 * flux_ix3() - sets BCs on L-x3 (left edge) of grid used in 3D
 * flux_ox3() - sets BCs on R-x3 (right edge) of grid used in 3D
 * grav_pot2() - gravitational potential for 2D problem (accn in Y)
 * grav_pot3() - gravitational potential for 3D problem (accn in Z)
 *============================================================================*/

static Real a,xi,Ftot,g,loz,delta,delad,Ttop,Tbot,Ptop,Pbot,T1,T0,P1,P0,ksmooth,minF;

void flux_ix2(GridS *pGrid);
void flux_ox2(GridS *pGrid);
void freeze_ix2(GridS *pGrid);
void freeze_ox2(GridS *pGrid);
void constT_ix2(GridS *pGrid);
void constT_ox2(GridS *pGrid);

void flux_ix3(GridS *pGrid);
void flux_ox3(GridS *pGrid);


Real heatcond_prof(const Real dens, const Real Press, const Real x1, const Real x2, const Real x3);
static Real grav_pot1(const Real x1, const Real x2, const Real x3);
static Real grav_pot2(const Real x1, const Real x2, const Real x3);
static Real grav_pot3(const Real x1, const Real x2, const Real x3);
static Real eix(Real x);
static Real e1xb(Real x); 
Real Tfunc(Real z);
Real Pfunc(Real z);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;
  long int iseed = -1;
  Real amp,x1,x2,x3,lx,ly,lz,rhoh,L_rot,fact;
#ifdef MHD
  Real b0,angle;
#endif
  Real Tval, Pval, Rhoval, Eval;
  int ixs, jxs, kxs;

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));


  kappa_iso = par_getd_def("problem","kappa_iso",1.0);
  a = par_getd_def("problem","a",3.5);
  xi = par_getd_def("problem","xi",3.0);
  Ftot = par_getd_def("problem","Ftot",1e-3);
  g = par_getd_def("problem","g",1.0);
  loz = par_getd_def("problem","loz",0.1);
  delta = par_getd_def("problem","delta",1e-3);
  ksmooth = par_getd_def("problem","ksmooth",0.05);
  minF = par_getd_def("problem","minF",0.1);

/* 2D PROBLEM --------------------------------------------------------------- */
/* Initialize two fluids with interface at y=0.0.  Pressure scaled to give a
 * sound speed of 1 at the interface in the light (lower, d=1) fluid 
 * Perturb V2 using single (iprob=1) or multiple (iprob=2) mode 
 */


  delad = 1.-1./Gamma;
  Ttop = xi/delad;
  T1 = Ttop + log((1+a)/(1+a*loz))/a;
  T0 = T1 + (1. + delta/delad)*(1+2*loz);

  Ptop = xi;
  P1 = Ptop *exp( (1+a)/delad *exp(a*Ttop)*(eix(-a*T1)-eix(-a*Ttop)));
  P0 = P1 *pow(T0/T1,1./(delad+delta));


  Tbot = Tfunc(-1.0);
  Pbot = Pfunc(-1.0);

  printf("Delad:%lg\nDelta:%lg\nTtop:%lg\nPtop:%lg\nT1:%lg\nT0:%lg\nP1:%lg\nP0:%lg\n",delad,delta,Ttop,Ptop,T1,T0,P1,P0);

//  FILE *f = fopen("/Users/zeus/Athena-Cversion/bin/ics.dat","w");

  if (pGrid->Nx[2] == 1) {
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      cc_pos(pGrid,1,j,k,&x1,&x2,&x3);
      Tval = Tfunc(x2);
      Pval = Pfunc(x2);
      Rhoval = Pval/(delad*Tval);
 //     fprintf(f,"%lg\t%lg\t%lg\t%lg\n",x2,Tval,Pval,Rhoval);
      for (i=is; i<=ie; i++) {
	    pGrid->U[k][j][i].d = Rhoval;
        pGrid->U[k][j][i].E = Pval/(Gamma-1);
	    pGrid->U[k][j][i].M1 = 0.0;
	    pGrid->U[k][j][i].M2 = 0.0;
	    pGrid->U[k][j][i].M3 = 0.0;

	    pGrid->U[k][j][i].E+=0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
	    pGrid->U[k][j][i].E+=0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
	    pGrid->U[k][j][i].E+=0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
	// Potential in E?
      }
    }
  }
//  fclose(f);

/* Enroll gravitational potential to give acceleration in y-direction for 2D
 * Use special boundary condition routines.  In 2D, gravity is in the
 * y-direction, so special boundary conditions needed for x2
*/

  StaticGravPot = grav_pot2;
 // CoolingFunc = AtmosCooling2D;
  HeatCondFunc = heatcond_prof;
 
  bvals_mhd_fun(pDomain, left_x2,  constT_ix2);
  bvals_mhd_fun(pDomain, right_x2, flux_ox2);

  } /* end of 2D initialization  */

/* 3D PROBLEM ----------------------------------------------------------------*/
/* Initialize two fluids with interface at z=0.0
 * Pressure scaled to give a sound speed of 1 at the interface
 * in the light (lower, d=1) fluid
 * iprob = 1 -- Perturb V3 using single mode
 * iprob = 2 -- Perturb V3 using multiple mode
 * iprob = 3 -- B in light fluid only, with multimode perturbation
 * iprob = 4 -- B rotated by "angle" at interface, multimode perturbation
 */

  if (pGrid->Nx[2] > 1) {
      printf("3D INIT!!\n");
  for (k=ks; k<=ke; k++) {
    cc_pos(pGrid,1,1,k,&x1,&x2,&x3);
      Tval = Tfunc(x3);
      Pval = Pfunc(x3);
      Rhoval = Pval/(delad*Tval);
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

        Rhoval = Pval/(delad*Tval);
        pGrid->U[k][j][i].d = Rhoval;
        pGrid->U[k][j][i].E = Pval/(Gamma-1);
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;

        pGrid->U[k][j][i].E+=0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].E+=0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
        pGrid->U[k][j][i].E+=0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
      }
    }
  }

/* Enroll gravitational potential to give accn in z-direction for 3D
 * Use special boundary condition routines.  In 3D, gravity is in the
 * z-direction, so special boundary conditions needed for x3
 */

  StaticGravPot = grav_pot3;
 // CoolingFunc = AtmosCooling3D;
  bvals_mhd_fun(pDomain, left_x3,  flux_ix3);
  bvals_mhd_fun(pDomain, right_x3, flux_ox3);

  } /* end of 3D initialization */

  return;
}
Real heatcond_prof(const Real dens, const Real Press, const Real x1, const Real x2, const Real x3) {
 
  Real zcval = (1.-minF)/a;


  Real result;
  Real logfac;

  Real xval = (x2 - zcval)/ksmooth;
  Real xmval = (-1. - zcval)/ksmooth;
  Real x1val = (1-x2 - zcval)/ksmooth;
  Real xm1val = (1.-(-1.) - zcval)/ksmooth;

  logfac = (1 + exp(xval))*(1+exp(-x1val));
  logfac /= (1+exp(xmval))*(1+exp(-xm1val));
  result =  1. - a*x2 + a*ksmooth*log(logfac);
 
  return result*Ftot;

}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll special boundary value functions,
 *    and initialize gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;

  if (pM->Nx[2] == 1) {
    StaticGravPot = grav_pot2;
    CoolingFunc = AtmosCooling2D;
    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x2,  constT_ix2);
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x2, flux_ox2);
      }
    }
  }
 
  if (pM->Nx[2] > 1) {
    StaticGravPot = grav_pot3; 
    CoolingFunc = AtmosCooling3D;
    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x3,  flux_ix3);
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x3, flux_ox3);
      }
    }
  }

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  \brief  Extracted from the Numerical Recipes in C (version 2) code.  
 *   Modified to use doubles instead of floats. - T. A. Gardiner - Aug. 12, 2003
 *   
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix2(GridS *pGrid)
 *  \brief constant flux, no-slip wall b.c 
 */

void freeze_ix2(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */
  Real Tval, Pval, Rhoval;
  Real x1,x2,x3,x10,x20,x30;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  cc_pos(pGrid,1,js,1,&x10,&x20,&x30);

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
        cc_pos(pGrid,1,js-j,k,&x1,&x2,&x3);
      Tval = Tfunc(x2);
      Pval = Pfunc(x2);
      Rhoval = Pval/(delad*Tval);
      for (i=il; i<=iu; i++) {

		pGrid->U[k][js-j][i].d = Rhoval; 
		pGrid->U[k][js-j][i].M1 = pGrid->U[k][js+(j-1)][i].M1;	// Zero gradient in horizontal velocities (stress free)	
		pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js+(j-1)][i].M2;	// Zero vertical velocity	
		pGrid->U[k][js-j][i].E = Pval/(Gamma-1) + (.5/Rhoval)*(SQR(pGrid->U[k][js-j][i].M1) + SQR(pGrid->U[k][js-j][i].M2)+SQR(pGrid->U[k][js-j][i].M3));
				
      }
    }
  }



  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x2 for 2D sims
 */

void freeze_ox2(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke, ku;
  int i,j,k,il,iu,jl,ju; /* i/j-lower/upper */
  Real Tval, Pval, Rhoval;
  Real x1,x2,x3,x10,x20,x30;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      cc_pos(pGrid,1,je+j,k,&x1,&x2,&x3);
      Tval = Tfunc(x2);
      Pval = Pfunc(x2);
      Rhoval = Pval/(delad*Tval);
      for (i=il; i<=iu; i++) {
		pGrid->U[k][je+j][i].d = Rhoval; 
		pGrid->U[k][je+j][i].M1 = pGrid->U[k][je-(j-1)][i].M1;	// Zero gradient in horizontal momenta (stress free)	
		pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je-(j-1)][i].M2;	// Zero vertical momentum 
		pGrid->U[k][je+j][i].E = Pval/(Gamma-1) + .5/Rhoval * (SQR(pGrid->U[k][je+j][i].M1) + SQR(pGrid->U[k][je+j][i].M2));
      }
    }
  }



  return;
}

void flux_ix2(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */
  Real Tval, Pval, Rhoval;
  Real Ts, Ps, Ds,Es,M1s,M2s;
  Real x1,x2,x3,x10,x20,x30;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  cc_pos(pGrid,1,js,1,&x10,&x20,&x30);

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
        cc_pos(pGrid,1,js-j,k,&x1,&x2,&x3);
      for (i=il; i<=iu; i++) {

		Ds = pGrid->U[k][js][i].d;				
		Es = pGrid->U[k][js][i].E;				
        M1s = pGrid->U[k][js][i].M1;
        M2s = pGrid->U[k][js][i].M2;
		Ps = (Gamma-1)*(Es - .5*(M1s*M1s+M2s*M2s)/Ds);
		Ts = Ps/(delad*Ds);

		Tval = Ts + log((1-a*x2)/(1-a*x20))/a; // Exact conductive T(z) 
		Pval = Ps*exp( (1-a*x20)/delad*exp(-a*Ts)*(eix(a*Tval)-eix(a*Ts))); // Exact hydrostatic pressure for T(z)
		Rhoval = Pval/(delad*Tval);

		pGrid->U[k][js-j][i].d = Rhoval; 
		pGrid->U[k][js-j][i].M1 = pGrid->U[k][js+(j-1)][i].M1;	// Zero gradient in horizontal velocities (stress free)	
		pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js+(j-1)][i].M2;	// Zero vertical velocity	
		pGrid->U[k][js-j][i].E = Pval/(Gamma-1) + (.5/Rhoval)*(SQR(pGrid->U[k][js-j][i].M1) + SQR(pGrid->U[k][js-j][i].M2)+SQR(pGrid->U[k][js-j][i].M3));
				
      }
    }
  }
  return;
}


void flux_ox2(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke, ku;
  int i,j,k,il,iu,jl,ju; /* i/j-lower/upper */
  Real Tval, Pval, Rhoval;
  Real Ts, Ps, Ds,Es;
  Real x1,x2,x3,x10,x20,x30;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  cc_pos(pGrid,1,je,1,&x10,&x20,&x30);
  x20 = 1-x20;
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
        cc_pos(pGrid,i,je+j,k,&x1,&x2,&x3);
		x2 = 1.-x2;
      for (i=il; i<=iu; i++) {
		Ds = pGrid->U[k][je][i].d;				
		Es = pGrid->U[k][je][i].E;				
		Ps = (Gamma-1)*(Es - .5*(SQR(pGrid->U[k][je][i].M1)+ SQR(pGrid->U[k][je][i].M2)+SQR(pGrid->U[k][je][i].M3))/Ds);
		Ts = Ps/(delad*Ds);

		Tval = Ts - log((1-a*x2)/(1-a*x20))/a; // Exact conductive T(z) 
		Pval = Ps*exp( (1-a*x20)/delad*exp(a*Ts)*(eix(-a*Tval)-eix(-a*Ts))); // Exact hydrostatic pressure for T(z)
		Rhoval = Pval/(delad*Tval);

		pGrid->U[k][je+j][i].d = Rhoval; 
		pGrid->U[k][je+j][i].M1 = pGrid->U[k][je-(j-1)][i].M1;	// Zero gradient in horizontal momenta (stress free)	
		pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je-(j-1)][i].M2;	// Zero vertical momentum 
		pGrid->U[k][je+j][i].E = Pval/(Gamma-1) + .5/Rhoval * (SQR(pGrid->U[k][je+j][i].M1) + SQR(pGrid->U[k][je+j][i].M2));
      }
    }
  }
  return;
}
void constT_ix2(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */
  Real Tval, Pval, Rhoval;
  Real Ts, Ps, Ds,Es,M1s,M2s;
  Real plaw;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {

		Ds = pGrid->U[k][js][i].d;				
		Es = pGrid->U[k][js][i].E;				
        M1s = pGrid->U[k][js][i].M1;
        M2s = pGrid->U[k][js][i].M2;
		Ps = (Gamma-1)*(Es - .5*(M1s*M1s+M2s*M2s)/Ds);
		Ts = Ps/(delad*Ds);

        
        Rhoval = pGrid->U[k][js+(j-1)][i].d;
        Pval = (Gamma-1)*(pGrid->U[k][js+(j-1)][i].E - .5*(SQR(pGrid->U[k][js+(j-1)][i].M1)+SQR(pGrid->U[k][js+(j-1)][i].M2))/Rhoval);

        Tval = 2*Tbot - Pval/(delad*Rhoval); // Linear profile of T through boundary

        plaw = -g/delad;
        plaw /= (Ts - Tbot)*2/pGrid->dx2;

        Pval = Ps * pow( Tval/Ts,plaw); // Hydrostatic Pressure given linear T profile
        Rhoval = Ds * pow( Tval/Ts,plaw-1); // Hydrostatic Density given EOS 

		pGrid->U[k][js-j][i].d = Rhoval; 
		pGrid->U[k][js-j][i].M1 = pGrid->U[k][js+(j-1)][i].M1;	// Zero gradient in horizontal velocities (stress free)	
		pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js+(j-1)][i].M2;	// Zero vertical velocity	
		pGrid->U[k][js-j][i].E = Pval/(Gamma-1) + (.5/Rhoval)*(SQR(pGrid->U[k][js-j][i].M1) + SQR(pGrid->U[k][js-j][i].M2)+SQR(pGrid->U[k][js-j][i].M3));
				
      }
    }
  }
  return;
}
void constT_ox2(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke, ku;
  int i,j,k,il,iu,jl,ju; /* i/j-lower/upper */
  Real Tval, Pval, Rhoval;
  Real Ts, Ps, Ds,Es,M1s,M2s;
  Real plaw;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
		Ds = pGrid->U[k][je][i].d;				
		Es = pGrid->U[k][je][i].E;				
        M1s = pGrid->U[k][je][i].M1;
        M2s = pGrid->U[k][je][i].M2;
		Ps = (Gamma-1)*(Es - .5*(M1s*M1s+M2s*M2s)/Ds);
		Ts = Ps/(delad*Ds);

        
        Rhoval = pGrid->U[k][je-(j-1)][i].d;
        Pval = (Gamma-1)*(pGrid->U[k][je-(j-1)][i].E - .5*(SQR(pGrid->U[k][je-(j-1)][i].M1)+SQR(pGrid->U[k][je-(j-1)][i].M2))/Rhoval);

        Tval = 2*Ttop - Pval/(delad*Rhoval); // Linear profile of T through boundary

        plaw = -g/delad;
        plaw /= (Ttop - Ts)*2/pGrid->dx2;

        Pval = Ps * pow( Tval/Ts,plaw); // Hydrostatic Pressure given linear T profile
        Rhoval = Ds * pow( Tval/Ts,plaw-1); // Hydrostatic Density given EOS 


		pGrid->U[k][je+j][i].d = Rhoval; 
		pGrid->U[k][je+j][i].M1 = pGrid->U[k][je-(j-1)][i].M1;	// Zero gradient in horizontal momenta (stress free)	
		pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je-(j-1)][i].M2;	// Zero vertical momentum 
		pGrid->U[k][je+j][i].E = Pval/(Gamma-1) + .5/Rhoval * (SQR(pGrid->U[k][je+j][i].M1) + SQR(pGrid->U[k][je+j][i].M2));
      }
    }
  }
  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix3(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x3 for 2D sims
 */

void flux_ix3(GridS *pGrid)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
  Real Tval, Pval, Rhoval;
  Real Ts, Ps, Ds,Es;
  Real x1,x2,x3,x10,x20,x30;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  cc_pos(pGrid,1,1,ks,&x10,&x20,&x30);
  for (k=1; k<=nghost; k++) {
    cc_pos(pGrid,1,1,k,&x1,&x2,&x3);
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
		Ds = pGrid->U[ks][j][i].d;				
		Es = pGrid->U[ks][j][i].E;				
		Ps = (Gamma-1)*(Es - .5*(SQR(pGrid->U[ks][j][i].M1) + SQR(pGrid->U[ks][j][i].M2)+SQR(pGrid->U[ks][j][i].M3))/Ds);
		Ts = Ps/(delad*Ds);

		Tval = Ts + log((1-a*x3)/(1-a*x30))/a; // Exact conductive T(z) 
		Pval = Ps*exp( (1-a*x30)/delad*exp(-a*Ts)*(eix(a*Tval)-eix(a*Ts))); // Exact hydrostatic pressure for T(z)
		Rhoval = Pval/(delad*Tval);

		pGrid->U[ks-k][j][i].d = Rhoval; 
		pGrid->U[ks-k][j][i].M1 = pGrid->U[ks+(k-1)][j][i].M1;	// Zero gradient in horizontal velocities (stress free)	
		pGrid->U[ks-k][j][i].M2 = pGrid->U[ks+(k-1)][j][i].M2;	// Zero gradient in horizontal velocities (stress free)		
		pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks+(k-1)][j][i].M3;	// Zero vertical velocity	
		pGrid->U[ks-k][j][i].E = Pval/(Gamma-1) + .5/Rhoval * (SQR(pGrid->U[ks-k][j][i].M1) + SQR(pGrid->U[ks-k][j][i].M2)+SQR(pGrid->U[ks-k][j][i].M3));

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x3 for 3D sims
 */

void flux_ox3(GridS *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k ,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
  Real Tval, Pval, Rhoval;
  Real Ts, Ps, Ds,Es;
  Real x1,x2,x3,x10,x20,x30;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }
  cc_pos(pGrid,1,1,ke,&x10,&x20,&x30);
  x30 = 1-x30;
  for (k=1; k<=nghost; k++) {
        cc_pos(pGrid,1,1,k,&x1,&x2,&x3);
		x3 = 1.-x3;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
		Ds = pGrid->U[ke][j][i].d;				
		Es = pGrid->U[ke][j][i].E;				
		Ps = (Gamma-1)*(Es - .5*(SQR(pGrid->U[ke][j][i].M1)+ SQR(pGrid->U[ke][j][i].M2)+SQR(pGrid->U[ke][j][i].M3))/Ds);
		Ts = Ps/(delad*Ds);

		Tval = Ts - log((1-a*x3)/(1-a*x30))/a; // Exact conductive T(z) 
		Pval = Ps*exp( (1-a*x30)/delad*exp(a*Ts)*(eix(-a*Tval)-eix(-a*Ts))); // Exact hydrostatic pressure for T(z)
		Rhoval = Pval/(delad*Tval);

		pGrid->U[ke+k][j][i].d = Rhoval; 
		pGrid->U[ke+k][j][i].M1 = pGrid->U[ke-(k-1)][j][i].M1;	// Zero gradient in horizontal momenta (stress free)	
		pGrid->U[ke+k][j][i].M2 = pGrid->U[ke-(k-1)][j][i].M2;	// Zero vertical momentum 
		pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke-(k-1)][j][i].M3;	// Zero vertical momentum 
		pGrid->U[ke+k][j][i].E = Pval/(Gamma-1) + .5/Rhoval * (SQR(pGrid->U[ke+k][j][i].M1) + SQR(pGrid->U[ke+k][j][i].M2)+SQR(pGrid->U[ke+k][j][i].M3));


      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot2(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 1.0
 */

static Real grav_pot1(const Real x1, const Real x2, const Real x3)
{
  return g*x1;
}
static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  return g*x2;
}
/*! \fn static Real grav_pot3(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */
static Real grav_pot3(const Real x1, const Real x2, const Real x3)
{
  return g*x3;
}


static Real e1xb(Real x) {
  int maxit;
  Real eps,euler,fpmin,fpmax;
  eps=6e-8;
  euler=.57721566;
  maxit=20;
  fpmin=1e-30;
  fpmax=1e30;
  int k;
  Real sum,term;

  if (x == 0.0) return fpmax;
  if (x <= 1.) {
    sum = 1.;
    term = 1.;
    k=1;
    while ((fabs(term)>fabs(sum)*eps)&&(k<=maxit)){
      term = -term*k*x/( (k+1.)*(k+1.));
      sum = sum + term;
      k++;
    }
    return -euler - log(x) + x*sum;
  }

  term = 0.;
  for (k=maxit+(int)floor(80/x);k>=1;k--) {
    term = k/(1. + k/(x+term));
  }
  return exp(-x) / (x+term);
}

static Real eix(Real x) {
  int maxit;
  Real eps,euler,fpmin,fpmax;
  eps=6e-8;
  euler=.57721566;
  maxit=100;
  fpmin=1e-30;
  fpmax=1e30;
  int k;
  Real sum,term;
  if (x == 0.0) return -fpmax;
  
  if (x < 0.0) return -e1xb(-x);
  if (fabs(x) <= 40.) {
    // Power series around x=0
  
    sum = 1.0;
    term =1.0;
    k = 1;
    while ((fabs(term/sum) > eps)&&(k<=maxit)) {
      term = term*k*x/( (k+1.)*(k+1.));
      sum += term;
      k++;
    }
    return euler + log(x) + x*sum;
  }
  
  // Asymptotic expansion (the series is not convergent)
  sum = 1.0;
  term = 1.0;
  
  for (k=1;k<=20;k++) {
    term = term*k/x;
    sum += term;
  }
  return exp(x)/x * sum;
}


Real Tfunc(Real z) {
    
    if (z < -loz) {
        return T0 + log((1-a*z)/(1+a*loz))/a;
    }

    if (z < 1+loz) {
        return T1 - (1. + delta/delad)*(z -(1+loz));
    }

    return Ttop - log((1-a*(1-z))/(1+a))/a;
}
Real Pfunc(Real z) {
    Real Tval = Tfunc(z);
    if (z < -loz) {
        return P0*exp((1+a*loz)/delad *exp(-a*T0)*(eix(a*Tval)-eix(a*T0)));
    }

    if (z < 1+loz) {
        return P1*pow(Tval/T1,1./(delad+delta));
    }

    return Ptop *exp((1+a)/delad *exp(a*Ttop)*(eix(-a*Tval)-eix(-a*Ttop)));
}
