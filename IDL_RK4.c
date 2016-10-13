# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
#include <pthread.h>
#include "idl_export.h"	
#include <unistd.h>

# define MAXERR 0.001 //Maximum relative error tolarance
# define MAXITER 10 //Maximum number of iterations
# define INITSTEPNUM 2 //initial number of steps

/*Thread function*/
void *solver_thread ( void *arg );

/* RK4 solver */
void rk4_solver ( double R_m[], double y[],  double t, double tout, double stepping);

/*Function calculating the derivative*/
void Deriv_func ( double t, double R_m[], double y[], double yp[], double stepping );

/* Adding arrays */
void array_add( double res[], double x1[], double x2[], int n_elements );
/* Multiplying array with constant  */
void array_multiply( double res[], double x[], double factor, int n_elements );
/* Finding maximum of the array  */
double array_max( double x[], int n_elements );
/* Getting absolute value of number */
double abs_val ( double x );

 
  int NLEV = 0; // Number of atomic levels
  int NPOINTS = 0; // Number of beam points
  int NBEAMS = 0; // Number of beams
  double *DIST; // Pointer to the dista matrix
  double *MATRIX; // Pointer to transition matrix
  double *STATES; // Pointer to the matrix of the populated states
  int nprocs=1; //Number of CPU-s
  
    

/******************************************************************************/

int gd_RK4( double *dista, double *C_m, int *n_levels, int *n_points, int *n_beams, double *pop, int *cpu_num )


/******************************************************************************/
/*
  Purpose: To calculate the evolution of atomic states along the given beams.

  Discussion: This calls the individual threads to solve the rate equations.

  Modified: 2011.07.14.

  Author: David Guszejnov

*/

{ 	

  int i;
  int *ind; // thread indices
  pthread_t *pth;	// thread identifier

/*Setting global variables from input*/
  NLEV = *n_levels;
  NPOINTS = *n_points;
  NBEAMS = *n_beams;
  DIST = dista;
  MATRIX = C_m;
  STATES = pop;
  nprocs =*cpu_num;

  
//    printf("CPU number %d", nprocs);
    
  /*No need for too many threads*/
  if(nprocs>3){
  if(NBEAMS < 60){nprocs=3;}
  if(NBEAMS < 40){nprocs=2;}
  if(NBEAMS < 20){nprocs=1;}
}

  /*So that others can use the system too*/  
  if(nprocs==*cpu_num){
    if(nprocs>2){
	 nprocs -= 1;
		 }
	 }
 //nprocs = 1;	//To test single threaded mode 
  
  //Pointer allocation for threads
  pth = ( pthread_t * ) malloc ( nprocs * sizeof ( pthread_t ) );
     
  //printf ( "\n Runge-Kutta integration started with %d thread(s) \n", nprocs );
  
  //Index array
  ind = ( int * ) malloc ( nprocs * sizeof ( int ) );
    
  /*Starting threads*/
  for ( i = 0; i < nprocs; i++ ){
  	  
   ind[i]=i;
  pthread_create(&pth[i],NULL,solver_thread,&ind[i]);
}

/*Waiting for threads*/
for ( i = 0; i < nprocs; i++ ){
 (void) pthread_join(pth[i], NULL);
}

 //printf( "\n Runge-Kutta integration finished \n");
 
//Cleanup
free(pth);
free(ind);

  return 1;
}

//*************************************************************************************

void *solver_thread ( void *arg){
  int i,j,threadID, start_beam, end_beam, n_step;
  double t,t_out,*y,*R_m,*D_m,*loc_pop, stepping;
  
  /*Thread specific data*/
  threadID=*((int*)arg);
   
  start_beam=NBEAMS/nprocs*threadID; //Starting beam index
  end_beam=NBEAMS/nprocs*(threadID+1)-1; //Ending beam index
   
  /*For last thread*/
  if((threadID+1)==nprocs){end_beam=NBEAMS-1;}

  
  //printf("Thread %d started for beams from %d to %d \n",threadID,start_beam, end_beam);

   /*Local arrays*/
  y = ( double * ) malloc ( NLEV * sizeof ( double ) );
  loc_pop = ( double * ) malloc ( (end_beam-start_beam+1)*NPOINTS*NLEV * sizeof ( double ) );
  R_m = ( double * ) malloc ( NPOINTS*NLEV*NLEV * sizeof ( double ) );
  D_m = ( double * ) malloc ( NPOINTS * sizeof ( double ) );

  /*For all beams*/
 for ( j = start_beam; j <= end_beam; j++ ){
  /*Set the transition matrix for the given beam*/
   for ( i = 0; i < NPOINTS*NLEV*NLEV; i++ ){
		   	   R_m[i]=MATRIX[i+j*NPOINTS*NLEV*NLEV];
			   }

  /*Set the distance array for the given beam*/
	for ( i = 0; i < NPOINTS; i++ ){
		D_m[i]=DIST[i+j*NPOINTS];
		}

  /*Defining stepping (assuming equidsitant grid)*/
  stepping=(D_m[NPOINTS-1]-D_m[0])/(double)(NPOINTS);
  
  /*Setting starting population*/
  y[0] = 1.0;
  STATES[NLEV*NPOINTS*j] = 1.0; //base pop
    for ( i = 1; i < NLEV; i++ ){
       y[i] = 0.0;
       STATES[NLEV*NPOINTS*j+i] = 0.0; //starting pop
     }
     

  /*Starting the integration*/
  for ( n_step = 1; n_step < NPOINTS ; n_step++ )
  {   
	t=D_m[n_step-1]; //Startpoint of calculation
	t_out=D_m[n_step];  //Endpoint of calculation
	/*Runge-Kutta sorlver started*/
    rk4_solver ( R_m, y, t, t_out, stepping);
    //printf("y0 %f at point %d\n",y[0],n_step);
    /*Storing results*/
     for ( i = 0; i < NLEV; i++ ){
       STATES[NLEV*NPOINTS*j+n_step*NLEV+i] = y[i];
	   }
  }
}

//Cleanup
free(y);
free(loc_pop);
free(R_m);
free(D_m);

  return;
}



void Deriv_func ( double t, double R_m[], double y[], double yp[], double stepping )

/******************************************************************************/
/*
  Purpose:

    Deriv_func evaluates the derivative for the ODE.


  Parameters:

    Input, double T, the value of the independent variable.

    Input, double Y(NLEV), the value of the dependent variable.
    
    Input, double R_m(NLEV x NPOINTS), the coefficients of the differential equation
    
    Input, double stepping, distance between two points

    Output double YP(NLEV), the value of the derivative dY(1:NEQN)/dT.
*/
{
  	int i,j,ind,step_num;
  	double fact;
	
   step_num=(int)(t/stepping); //Getting the stepnumber
   fact=t/stepping-(double)(step_num); //Factor for linear interpolation
   //Correcting for the ending
   if(step_num>=(NPOINTS-1)){
   						 step_num=NPOINTS-2;
   						 fact=1.0;
							} 
   

   /*Calculating derivative*/ 
  for(i=0;i<NLEV;i++){
          yp[i]=0;  
  		  for(j=0;j<NLEV;j++){
 			 ind=step_num*NLEV*NLEV+i*NLEV+j; //Index of the array
  				yp[i]+=((1-fact)*R_m[ind]+fact*R_m[ind+NLEV*NLEV])*y[j]; 
						}				  
							}
								 					  
  
  return;
}


/******************************************************************************/

void rk4_solver ( double R_m[], double y[], double t, double tout, double stepping )


/******************************************************************************/
/*
  Purpose:

    Carries out the 4th order Runge-Kutta integration between t and tout using error tolerance of MAXERR
*/
{

  int i,j; //loop integer
  int iteration=0; //iteration number
  int nsteps; //number of integration steps
  double timestep; //time step
  double tcurr; //current time
  double *yt, *yold, *temp, *k1, *k2, *k3, *k4; //intermediate results
  double *relerr, maxrelerr=1.0; //Relative error of calculation
  //Allocating memory for arrays
  yt = ( double * ) malloc ( NLEV * sizeof ( double ) );
  yold = ( double * ) malloc ( NLEV * sizeof ( double ) );
  temp = ( double * ) malloc ( NLEV * sizeof ( double ) );
  k1 = ( double * ) malloc ( NLEV * sizeof ( double ) );
  k2 = ( double * ) malloc ( NLEV * sizeof ( double ) );
  k3 = ( double * ) malloc ( NLEV * sizeof ( double ) );
  k4 = ( double * ) malloc ( NLEV * sizeof ( double ) );
  relerr = ( double * ) malloc ( NLEV * sizeof ( double ) );
  
  //Initial values
  nsteps=INITSTEPNUM;
  
  //Starting loop
  while ((maxrelerr > MAXERR) && (iteration < MAXITER)) {
  		  nsteps*=2;
  		  tcurr=t; //resetting timer
  		  timestep=(tout-t)/(double)nsteps;
  		  iteration+=1;
  		 // printf("y: %f %f %f %f \n",y[0],y[1],y[2],y[3]);
  		  //Resetting y vector
  		    for(i=0;i<NLEV;i++){
					yt[i]=y[i];
					}
  		  //printf("RK4 starting Iteration %d \n",iteration);
  		for(i=0;i<nsteps;i++){
							  Deriv_func ( tcurr , R_m, yt, k1, stepping ); //calculating k1
							  array_multiply( k1, k1, timestep, NLEV );
							 
							  	array_multiply( temp, k1, 0.5, NLEV ); // getting 0.5*k1
							  	array_add( temp, yt, temp, NLEV ); // getting y+0.5*k1
							  Deriv_func ( (tcurr+0.5*timestep) , R_m, temp, k2, stepping ); //calculating k2
							  array_multiply( k2, k2, timestep, NLEV ); 
							  	array_multiply( temp, k2, 0.5, NLEV ); // getting 0.5*k2
							  	array_add( temp, yt, temp, NLEV ); // getting y+0.5*k2
							  Deriv_func ( (tcurr+0.5*timestep) , R_m, temp, k3, stepping ); //calculating k3
							  array_multiply( k3, k3, timestep, NLEV ); 
							  	array_add( temp, yt, k3, NLEV ); // getting y+k3
							  Deriv_func ( (tcurr+timestep) , R_m, temp, k4, stepping ); //calculating k4
							  array_multiply( k4, k4, timestep, NLEV );
							  //Getting next y values
	   						  for(j=0;j<NLEV;j++){
							  					  yt[j]+=(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6.0;
												  }		  						
							  tcurr+=timestep; //advancing time
							  }
		  //Calculating relative error
		  if (iteration > 1){
		  	 			for(i=0;i<NLEV;i++){
											relerr[i]=abs_val(yold[i]/yt[i]-1.0);
											}
						maxrelerr= array_max( relerr, NLEV );
						}
			//printf("Relerr calculated in iteration %d \n",iteration);
		  //Storing values of iteration
		  for(i=0;i<NLEV;i++){
							yold[i]=yt[i];
							}
						//ending while segment	
							}
//Storing results
for(i=0;i<NLEV;i++){
					y[i]=yt[i];
					}	

//Cleanup
free(yt);
free(yold);
free(temp);
free(k1);
free(k2);
free(k3);
free(k4);
free(relerr);

return;           

      }
          
 	



void array_add( double res[], double x1[], double x2[], int n_elements )
{
 	 //adds up two arrays
 	   int i;	     
 	   for(i=0;i<n_elements;i++){
	   							 res[i]=x1[i]+x2[i];	   
 	   
	   }
	   return;
	   }
void array_multiply( double res[], double x[], double factor, int n_elements )
{
 	 //multiplies array with constant
 	   int i;	     
 	   for(i=0;i<n_elements;i++){
	   							 res[i]=x[i]*factor;	   
 	   
	   }
	   return;
	   }
	   
double array_max( double x[], int n_elements )
	   {
 	 //Finds maximum of the array
 	   int i;
	   double max=x[0];	     
 	   for(i=0;i<n_elements;i++){
			  if ( max < x[i] )	{
			  	   max=x[i];
			  	 }		   	   
	   }
	   return max;
	   }
	   
	   
double abs_val ( double x )
//Gets absolute value of number
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
