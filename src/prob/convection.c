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
 * reflect_ix2() - sets BCs on L-x2 (left edge) of grid used in 2D
 * reflect_ox2() - sets BCs on R-x2 (right edge) of grid used in 2D
 * reflect_ix3() - sets BCs on L-x3 (left edge) of grid used in 3D
 * reflect_ox3() - sets BCs on R-x3 (right edge) of grid used in 3D
 * grav_pot2() - gravitational potential for 2D problem (accn in Y)
 * grav_pot3() - gravitational potential for 3D problem (accn in Z)
 *============================================================================*/

static double ran2(long int *idum);
static void reflect_ix2(GridS *pGrid);
static void reflect_ox2(GridS *pGrid);
static void reflect_ix3(GridS *pGrid);
static void reflect_ox3(GridS *pGrid);
static Real grav_pot2(const Real x1, const Real x2, const Real x3);
static Real grav_pot3(const Real x1, const Real x2, const Real x3);
static Real eix(Real x);
static Real e1xb(Real x); 

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
  Real T0,T1,P0,P1,Ttop,Ptop;
  Real Ftot,xi,dlnkdz,g,delad;
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

/* Read perturbation amplitude, problem switch, background density */
  amp = par_getd("problem","amp");
  iprob = par_geti("problem","iprob");
  rhoh  = par_getd_def("problem","rhoh",3.0);
/* Distance over which field is rotated */
  L_rot  = par_getd_def("problem","L_rot",0.0);

/* Read magnetic field strength, angle [should be in degrees, 0 is along +ve
 * X-axis (no rotation)] */
#ifdef MHD
  b0 = par_getd("problem","b0");
  angle = par_getd("problem","angle");
  angle = (angle/180.)*PI;
#endif

  dlnkdz = par_getd("problem","dlnkdz");
  xi = par_getd("problem","xi");
  Ftot = par_getd("problem","Ftot");
  g = par_getd("problem","g");
  loz = par_getd("problem","loz");
  delta = par_getd("problem","delta");

/* 2D PROBLEM --------------------------------------------------------------- */
/* Initialize two fluids with interface at y=0.0.  Pressure scaled to give a
 * sound speed of 1 at the interface in the light (lower, d=1) fluid 
 * Perturb V2 using single (iprob=1) or multiple (iprob=2) mode 
 */


  delad = 1.-1./Gamma;
  Ttop = xi/delad;
  Ptop = xi;
  T1 = Ttop - log( (1+dlnkdz*loz)/(1+dlnkdz))/dlnkdz; 
  T0 = T1 + (delad+delta)*(1+2*loz)/delad;
  P1 = Ptop*exp( (1+dlnkdz)/delad * exp(dlnkdz*Ttop) *( eix(-dlnkdz*T1)-eix(-dlnkdz*Ttop))); 
  P0 = P1*pow(T0/T1,1./(delad+delta));

  if (pGrid->Nx[1] == 1) {
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	if (x1 < -loz) {
	  Tval = T0 + log( (1-dlnkdz*x1)/(1+dlnkdz*loz))/dlnkdz;
	  Pval = P0*exp( (1+dlnkdz*loz)/delad * exp(-dlnkdz*T0) *( eix(dlnkdz*Tval)-eix(dlnkdz*T0))); 

	}
	else if (x1 < 1+loz) {
	  Tval = T1 - g*(1 + delta/delad)*(x1 - (1+loz));
	  Pval = P1*pow(Tval/T1,1./(delad+delta));

	}
	else {
	  Tval = Ttop - log( (1-dlnkdz*(1-x1))/(1+dlnkdz))/dlnkdz;
	  Pval = Ptop*exp( (1+dlnkdz)/delad * exp(dlnkdz*Ttop) *( eix(-dlnkdz*Tval)-eix(-dlnkdz*Ttop))); 


	}

	Dval = Pval/(delad*Tval);

	Eval = Pval/(Gamma-1);



	pGrid->U[k][j][i].d = Dval;
        pGrid->U[k][j][i].E = Eval;
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

/* Enroll gravitational potential to give acceleration in y-direction for 2D
 * Use special boundary condition routines.  In 2D, gravity is in the
 * y-direction, so special boundary conditions needed for x2
*/

  StaticGravPot = grav_pot1;
  CoolingFunc = AtmosCool1D;
  if (pDomain->Disp[1] == 0) bvals_mhd_fun(pDomain, left_x1,  reflect_ix1);
  if (pDomain->MaxX[1] == pDomain->RootMaxX[1])
    bvals_mhd_fun(pDomain, right_x1, reflect_ox1);

  } /* end of 1D initialization  */
  if (pGrid->Nx[2] == 1) {
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	if (x2 < -loz) {
	  Tval = T0 + log( (1-dlnkdz*x2)/(1+dlnkdz*loz))/dlnkdz;
	  Pval = P0*exp( (1+dlnkdz*loz)/delad * exp(-dlnkdz*T0) *( expi(dlnkdz*Tval)-expi(dlnkdz*T0))); 

	}
	else if (x2 < 1+loz) {
	  Tval = T1 - g*(1 + delta/delad)*(x2 - (1+loz));
	  Pval = P1*pow(Tval/T1,1./(delad+delta));

	}
	else {
	  Tval = Ttop - log( (1-dlnkdz*(1-x2))/(1+dlnkdz))/dlnkdz;
	  Pval = Ptop*exp( (1+dlnkdz)/delad * exp(dlnkdz*Ttop) *( expi(-dlnkdz*Tval)-expi(-dlnkdz*Ttop))); 


	}

	Dval = Pval/(delad*Tval);

	Eval = Pval/(Gamma-1);



	pGrid->U[k][j][i].d = Dval;
        pGrid->U[k][j][i].E = Eval;
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

/* Enroll gravitational potential to give acceleration in y-direction for 2D
 * Use special boundary condition routines.  In 2D, gravity is in the
 * y-direction, so special boundary conditions needed for x2
*/

  StaticGravPot = grav_pot2;
  CoolingFunc = AtmosCool2D;
  if (pDomain->Disp[1] == 0) bvals_mhd_fun(pDomain, left_x2,  reflect_ix2);
  if (pDomain->MaxX[1] == pDomain->RootMaxX[1])
    bvals_mhd_fun(pDomain, right_x2, reflect_ox2);

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
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	if (x3 < -loz) {
	  Tval = T0 + log( (1-dlnkdz*x3)/(1+dlnkdz*loz))/dlnkdz;
	  Pval = P0*exp( (1+dlnkdz*loz)/delad * exp(-dlnkdz*T0) *( expi(dlnkdz*Tval)-expi(dlnkdz*T0))); 

	}
	else if (x3 < 1+loz) {
	  Tval = T1 - g*(1 + delta/delad)*(x3 - (1+loz));
	  Pval = P1*pow(Tval/T1,1./(delad+delta));

	}
	else {
	  Tval = Ttop - log( (1-dlnkdz*(1-x3))/(1+dlnkdz))/dlnkdz;
	  Pval = Ptop*exp( (1+dlnkdz)/delad * exp(dlnkdz*Ttop) *( expi(-dlnkdz*Tval)-expi(-dlnkdz*Ttop))); 


	}

	Dval = Pval/(delad*Tval);

	Eval = Pval/(Gamma-1);



	pGrid->U[k][j][i].d = Dval;
        pGrid->U[k][j][i].E = Eval;
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
  CoolingFunc = AtmosCool3D;

  if (pDomain->Disp[2] == 0) bvals_mhd_fun(pDomain, left_x3,  reflect_ix3);
  if (pDomain->MaxX[2] == pDomain->RootMaxX[2])
    bvals_mhd_fun(pDomain, right_x3, reflect_ox3);

  } /* end of 3D initialization */

  return;
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
    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x2,  reflect_ix2);
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x2, reflect_ox2);
      }
    }
  }
 
  if (pM->Nx[2] > 1) {
    StaticGravPot = grav_pot3;
    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x3,  reflect_ix3);
        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x3, reflect_ox3);
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

static void flux_ix2(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */
  Real Tval, Pval, Rhoval
  Real delad = 1. - 1./Gamma;
	Real x1,x2,x3,x10,x20,x30;
  Real dlnkdz = 3.5;

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  cc_pos(pGrid,1,js,1,&x10,&x20,&x30);

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

				D0 = pGrid->U[k][js][i].d;				
				E0 = pGrid->U[k][js][i].E;				
				P0 = (Gamma-1)*(E0 - .5*(SQR(pGrid->U[k][js][i].M1) SQR(pGrid->U[k][js][i].M2)-SQR(pGrid->U[k][js][i].M3))/D0);
				T0 = P0/(delad*D0);

				Tval = T0 + log((1-dlnkdz*x2)/(1-dlnkdz*x20))/dlnkdz; // Exact conductive T(z) 
				Pval = P0*exp( (1-dlnkdz*x20)/delad*exp(-dlnkdz*T0)*(eix(dlnkdz*Tval)-eix(dlnkdz*T0))); // Exact hydrostatic pressure for T(z)
				Rhoval = Pval/(delad*Tval);

				pGrid->U[k][js-j][i].d = Rhoval; 
				pGrid->U[k][js-j][i].M1 = pGrid->U[k][js+(j-1)][i].M1;	// Zero gradient in horizontal velocities (stress free)	
				pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js+(j-1)][i].M2;	// Zero vertical velocity	
				pGrid->U[k][js-j][i].d = Pval/(Gamma-1) + .5/Rhoval * (SQR(pGrid->U[k][js-j][i].M1) + SQR(pGrid->U[k][js-j][i].M2)+SQR(pGrid->U[k][js-j][i].M3));
				
      }
    }
  }



  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x2 for 2D sims
 */

static void reflect_ox2(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke, ku;
  int i,j,k,il,iu,jl,ju; /* i/j-lower/upper */
  Real Tval, Pval, Rhoval
  Real delad = 1. - 1./Gamma;
	Real x1,x2,x3,x10,x20,x30;
  Real dlnkdz = 3.5;

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
      for (i=il; i<=iu; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
				x2 = 1.-x2;

				D0 = pGrid->U[k][je][i].d;				
				E0 = pGrid->U[k][je][i].E;				
				P0 = (Gamma-1)*(E0 - .5*(SQR(pGrid->U[k][je][i].M1) SQR(pGrid->U[k][je][i].M2)-SQR(pGrid->U[k][je][i].M3))/D0);
				T0 = P0/(delad*D0);

				Tval = T0 - log((1-dlnkdz*x2)/(1-dlnkdz*x20))/dlnkdz; // Exact conductive T(z) 
				Pval = P0*exp( (1-dlnkdz*x20)/delad*exp(dlnkdz*T0)*(eix(-dlnkdz*Tval)-eix(-dlnkdz*T0))); // Exact hydrostatic pressure for T(z)
				Rhoval = Pval/(delad*Tval);

				pGrid->U[k][je+j][i].d = Rhoval; 
				pGrid->U[k][je+j][i].M1 = pGrid->U[k][je-(j-1)][i].M1;	// Zero gradient in horizontal momenta (stress free)	
				pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je-(j-1)][i].M2;	// Zero vertical momentum 
				pGrid->U[k][je+j][i].d = Pval/(Gamma-1) + .5/Rhoval * (SQR(pGrid->U[k][je+j][i].M1) + SQR(pGrid->U[k][je+j][i].M2)+SQR(pGrid->U[k][je+j][i].M3));
      }
    }
  }



  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix3(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x3 for 2D sims
 */

static void reflect_ix3(GridS *pGrid)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ks-k][j][i]    =  pGrid->U[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks-k][j][i].M3; /* reflect 3-mom. */
        pGrid->U[ks-k][j][i].E +=
          pGrid->U[ks+(k-1)][j][i].d*0.1*(2*k-1)*pGrid->dx3/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks+(k-1)][j][i];
      }
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3(GridS *pGrid)
 *  \brief Special reflecting boundary functions in x3 for 3D sims
 */

static void reflect_ox3(GridS *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k ,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ke+k][j][i]    =  pGrid->U[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke+k][j][i].M3; /* reflect 3-mom. */
        pGrid->U[ke+k][j][i].E -=
          pGrid->U[ke-(k-1)][j][i].d*0.1*(2*k-1)*pGrid->dx3/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke-(k-1)][j][i];
      }
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot2(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 1.0
 */

static Real grav_pot1(const Real x1, const Real x2, const Real x3)
{
  return x1;
}
static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  return x2;
}
/*! \fn static Real grav_pot3(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */
static Real grav_pot3(const Real x1, const Real x2, const Real x3)
{
  return x3;
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
  
  if (x < 0.0) return -de1xb(-x);
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
