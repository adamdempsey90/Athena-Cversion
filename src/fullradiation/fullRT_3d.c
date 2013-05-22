#include "../copyright.h"
/*==============================================================================
 * FILE: formal_solution.c
 *
 * PURPOSE: integrator for the full radiation transfer in 2D 
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   formal_solution()  - interate formal solution until convergence
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef FULL_RADIATION_TRANSFER


static Real ******Divi = NULL;	/* temporary array to store flux for each array */
static Real ****flux = NULL;
static Real *tempimu = NULL;




/* calculate the upwind and downwind intensity for k,j,i, frequency ifr, octant l and angle n */
/* This function is only used to calculate for one ray */
/* in 3D, the octant is numbered as *
 *           | 
 *	1    |    0 
 * ---------------------
 *	3    |    2
 *	     |
 *----------------------
 *---------------------
 *           | 
 *	5    |    4 
 * ---------------------
 *	7    |    6
 *	     |
 *----------------------
 ********************/



void fullRT_3d(DomainS *pD)
{

	RadGridS *pRG=(pD->RadGrid);
	Real dt = (pD->Grid)->dt;
	int i, is, ie;
	int j, js, je;
	int k, ks, ke;
	is = pRG->is; ie = pRG->ie;
	js = pRG->js; je = pRG->je;
	ks = pRG->ks; ke = pRG->ke;

	int l, n, ifr, m;

	Real dx1, dx2, dx3, ds;
	dx1 = pRG->dx1;
	dx2 = pRG->dx2;
	dx3 = pRG->dx3;

	Real imu[5];
	Real vel;




	/* first, update the source functions */
	hydro_to_fullrad(pD);


	
	
	for(ifr=0; ifr<pRG->nf; ifr++){
		for(l=0; l<pRG->noct; l++){
			for(n=0; n<pRG->nang; n++){

				/* Now calculate the x flux */
				ds = dx1;
				for(k=ks; k<=ke; k++){	
					for(j=js; j<=je; j++){					
						if((l == 0) || (l == 2) || (l == 4) || (l == 6)){
							for(i=is; i<=ie+1; i++){
								/* This is only true for constant angle */
								/* Need to consider when advection with non-uniform velocity */
								vel = Crat;
								/* From small i to large i */
#ifdef SECOND_ORDER_PRIM
								for(m=0; m<3; m++)
									imu[m] = pRG->imu[ifr][l][n][k][j][i-2+m] * pRG->mu[l][n][k][j][i-2+m][0];

							

								flux_PLM(dt, ds, vel, imu, &(flux[ifr][l][n][i]));	
#else							/*------------------------------------------------------------------*/
								for(m=0; m<5; m++)
									imu[m] = pRG->imu[ifr][l][n][k][j][i-3+m] * pRG->mu[l][n][k][j][i-3+m][0];

								flux_PPM(dt, ds, vel, imu, &(flux[ifr][l][n][i]));
#endif							
													

							}/* end i */		
						}/* end l==0, 2, 4, 6*/
						else{
							for(i=is-1; i<=ie; i++){
								/* From large i to small i */
								vel = Crat;
#ifdef SECOND_ORDER_PRIM
								for(m=2; m>=0; m--)
									imu[m] = pRG->imu[ifr][l][n][k][j][i+2-m] * pRG->mu[l][n][k][j][i+2-m][0];

								flux_PLM(dt, ds, vel, imu, &(flux[ifr][l][n][i+1]));		

#else							/*------------------------------------------------------------------*/

								for(m=4; m>=0; m--)
									imu[m] = pRG->imu[ifr][l][n][k][j][i+3-m] * pRG->mu[l][n][k][j][i+3-m][0];

								flux_PPM(dt, ds, vel, imu, &(flux[ifr][l][n][i+1]));
#endif	

							} /* end i */
						}/* end l == 1, 3, 5, 7 */
	
		/* Now save the flux difference */
				
			   			for(i=is; i<=ie; i++){
							Divi[ifr][l][n][k][j][i] =  Crat * dt * (flux[ifr][l][n][i+1] - flux[ifr][l][n][i]) /(dx1);	
		   				}/* end i */
					}/* End j */
				}/* End k */

				/*---------------------------------------------------------------------*/

				ds = dx2;
			/* Now calculate the flux along j direction */
				for(k=ks; k<=ke; k++){
					for(i=is; i<=ie; i++){
						/* first save the data */
						for (j=0; j<=je+Radghost; j++){
							tempimu[j] = pRG->imu[ifr][l][n][k][j][i] * pRG->mu[l][n][k][j][i][1];
						}
					
						if((l == 0) || (l == 1) || (l == 4) || (l == 5)){
							for(j=js; j<=je+1; j++){
								vel = Crat;
						
#ifdef SECOND_ORDER_PRIM
								for(m=0; m<3; m++)
									imu[m] = tempimu[j-2+m];
							
								flux_PLM(dt, ds, vel, imu, &(flux[ifr][l][n][j]));	
#else
								for(m=0; m<5; m++)
									imu[m] = tempimu[j-3+m];

								flux_PPM(dt, ds, vel, imu, &(flux[ifr][l][n][j]));
#endif							

							}/* end j */
						}/* end l ==0, 1, 4, 5, */
						else{
							for(j=js-1; j<=je; j++){
								vel = Crat;
#ifdef SECOND_ORDER_PRIM
								for(m=2; m>=0; m--)
									imu[m] = tempimu[j+2-m];
			
								flux_PLM(dt, ds, vel, imu, &(flux[ifr][l][n][j+1]));	
#else
								for(m=4; m>=0; m--)
									imu[m] = tempimu[j+3-m];

								flux_PPM(dt, ds, vel, imu, &(flux[ifr][l][n][j+1]));
#endif				
							} /* end j */				
						}/* end l == 2, 3, 6, 7 */				
						
		/* Now save the flux difference */
		
		
		   				for(j=js; j<=je; j++){
							Divi[ifr][l][n][k][j][i] += (Crat * dt * (flux[ifr][l][n][j+1] - flux[ifr][l][n][j]) /(dx2));
						} /* end j */
					} /* Finish i */
				}/* Finish k */


				/*---------------------------------------------------------------------*/

				/* Now calculate the flux along x3 direction */
				ds = dx3;

				for(j=js; j<=je; j++){
					for(i=is; i<=ie; i++){
						/* first save the data */
						for (k=0; k<=ke+Radghost; k++){
							tempimu[k] = pRG->imu[ifr][l][n][k][j][i] * pRG->mu[l][n][k][j][i][2];
						}
					
						if((l == 0) || (l == 1) || (l == 2) || (l == 3)){
							for(k=ks; k<=ke+1; k++){
								vel = Crat;
						
#ifdef SECOND_ORDER_PRIM
								for(m=0; m<3; m++)
									imu[m] = tempimu[k-2+m];
							
								flux_PLM(dt, ds, vel, imu, &(flux[ifr][l][n][k]));	
#else
								for(m=0; m<5; m++)
									imu[m] = tempimu[k-3+m];

								flux_PPM(dt, ds, vel, imu, &(flux[ifr][l][n][k]));
#endif							

							}/* end k */
						}/* end l ==0, 1, 2, 3, */
						else{
							for(k=ks-1; k<=ke; k++){
								vel = Crat;
#ifdef SECOND_ORDER_PRIM
								for(m=2; m>=0; m--)
									imu[m] = tempimu[k+2-m];
			
								flux_PLM(dt, ds, vel, imu, &(flux[ifr][l][n][k+1]));	
#else
								for(m=4; m>=0; m--)
									imu[m] = tempimu[k+3-m];

								flux_PPM(dt, ds, vel, imu, &(flux[ifr][l][n][k+1]));
#endif				
							} /* end k */				
						}/* end l == 4, 5, 6, 7 */				
						
		/* Now save the flux difference */
		
		
		   				for(k=ks; k<=ke; k++){
							Divi[ifr][l][n][k][j][i] += (Crat * dt * (flux[ifr][l][n][k+1] - flux[ifr][l][n][k]) /(dx3));
						} /* end k */
					} /* Finish i */
				}/* Finish j */



			}/* end n */
		}/* end l */
	}/* end ifr */


	/* Now we have flux and source term for each array, update them together */
	for(ifr=0; ifr<pRG->nf; ifr++){
		for(l=0; l<pRG->noct; l++){
			for(n=0; n<pRG->nang; n++){
				for(k=ks; k<=ke; k++){
					for(j=js; j<=je; j++){
						for(i=is; i<=ie; i++){
							pRG->imu[ifr][l][n][k][j][i] = (pRG->heatcool[ifr][l][n][k][j][i] +  (pRG->imu[ifr][l][n][k][j][i]- Divi[ifr][l][n][k][j][i] + dt * Crat * pRG->R[ifr][k][j][i].Sigma[2] * pRG->R[ifr][k][j][i].J)/(1.0 + dt * Crat * pRG->R[ifr][k][j][i].Sigma[2]));
						}/* end i */
					}/* end j */
				}/* end k */
			}/* end n */
		}/* end l */
	}/* end ifr */
		

	/* Moments are updated in the main loop */
	

  return;
}



void fullRT_3d_destruct(void)
{

 

  if (Divi != NULL) free_6d_array(Divi);	
  if (flux != NULL) free_4d_array(flux);
  if(tempimu != NULL) free_1d_array(tempimu);
  return;
}


void fullRT_3d_init(RadGridS *pRG)
{

	
	int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1], nx3 = pRG->Nx[2];
	int nfr = pRG->nf, noct = pRG->noct, nang = pRG->nang;
	int nmax, nmaxjk;


	nmax = MAX(nx1,nx2);
	nmax = MAX(nmax,nx3);

	nmaxjk = MAX(nx2, nx3);

	if ((flux = (Real ****)calloc_4d_array(nfr, noct, nang, nmax+2*Radghost, sizeof(Real))) == NULL)
    		goto on_error;


	if ((Divi = (Real ******)calloc_6d_array(nfr, noct, nang, nx3+2*Radghost, nx2+2*Radghost, nx1+2*Radghost, sizeof(Real))) == NULL)
    		goto on_error;

	if ((tempimu = (Real *)calloc_1d_array(nmaxjk+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;

	return;

	on_error:
  	fullRT_3d_destruct();
  	ath_error("[fullRT_3d_init]: Error allocating memory\n");
  	return;

}
#endif /* FULL_RADIATION_TRANSFER */
