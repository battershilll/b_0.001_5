
//Two phase with diffusive interface

#include "vof.h"
#include "tracer.h"
#include "diffusion.h"

//#define G 9.81

const face vector av[];
scalar f3[];
scalar * tracers = {f3};
scalar f[];
scalar * interfaces  = {f}; // Tracer for the vol fraction
mgstats mgT;

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1. , mu3 = 0.;  // we initlise the 

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
	  alpha = alphav;
	    rho = rhov;
  if (mu1 || mu2)
	      mu = new face vector;
}


#ifndef rho
# define rho(f,f3) (clamp(f,0.,1.)*(clamp(f3,0.,1.)*(rho3 - rho1) + rho1 -rho2) + rho2) //void fraction (density) is a function of the volume fraction. can add on a third phase here. 
#endif
#ifndef mu
# define mu(f,f3)  (clamp(f,0.,1.)*(clamp(f3,0.,1.)*(mu3 - mu1) + mu1 - mu2) + mu2) //viscoisty is a function of the volumne fraction
#endif

		/*
//Acceleration event
event acceleration (i++) {
	  foreach_face(y)
		      av.y[] -= 9.81;
	    boundary ((scalar *){av});
}
*/

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event properties (i++) {
	//When using smearing of the density jump, we initialise sf with the vertex-average of f.

#ifndef sf
#if dimension <= 2
		  foreach()
		      sf[] = (4.*f[] + 
				      	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
					    	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
	  foreach()
		      sf[] = (8.*f[] +
				      	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
					    	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] + 
							    		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
											f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
						    	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
							    	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif 

#if TREE
	    sf.prolongation = refine_bilinear;
	      boundary ({sf});
#endif

	        
	        foreach_face() {
			    double ff = (f[] + f[-1])/2.;
			        alphav.x[] = fm.x[]/rho(ff,f3[]);
				    if (mu1 || mu2) {
					          face vector muv = mu;
						        muv.x[] = fm.x[]*mu(ff,f3[]);
							    }
				      }
		  foreach()
			      rhov[] = cm[]*rho(f[],f3[]);

#if TREE  
		    sf.prolongation = fraction_refine;
		      boundary ({sf});
#endif

}

//event for the diffusion of the concentration

event tracer_diffusion (i++) {
//#if dimension == 2
	//We want the diffusion vector to be dependent on the phase we are in

       const face vector D[] = {0.0001,0.0001};	
	diffusion (f3, dt, D);
	 boundary ({f3});


}




