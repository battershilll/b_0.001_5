#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
scalar f3[];
double rho3 = 1400; //denstity for dense salt water solution
double mu3 = 0.001;

//scalar * tracers = {f3};

//We “overload” the default by defining the mu() macro before including the code for two phases.

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(clamp(f3[],0.,1.)*(rho3 - rho1) + rho1 -rho2) + rho2) //void fraction (density) is a function of the volume fraction. can add on a third phase here
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(clamp(f3[],0.,1.)*(mu3 - mu1) + mu1 - mu2) + mu2) //viscoisty is a function of the volumne fraction
#endif

/*
#ifndef diff
# define diff(f3) (clamp(f3,0.,1)*0.1)
#endif*/

//#include "advection_Q.h"
#include "two-phase.h"
#include "mixtures.h"
//#include "tracer.h"
//#include "elementary_body.h"

#include "utils.h"


FILE * fp0;

/**
  The domain is eight units long, centered vertically. */

int main() {
	  L0 = 8.;
	    origin (0, 0);
	    rho1 = 1010; //water
	    rho2 = 100; //air 
	    mu1 = 0.001;
	    mu2 = 0.00001;
	    //mu3 = 0.0001;
	    init_grid(1024); //Level 10
	    //f.sigma = 73e-3;
	    //G.x = 0.;
	    //G.y = -9.81;
	    fp0=fopen("tracers_cons_final.dat", "w");
 	    const face vector av[] = {9.81*sin(0.268), -9.81*cos(0.268), 0};
	    a = av;
	    DT = 0.001;
	   // TOLERANCE = 1e-;
	    run(); 
}


u.n[top] = dirichlet (0);
u.t[top] = dirichlet (0);
u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
u.n[left] = dirichlet (0);
u.t[left] = dirichlet (0);
u.n[right] = dirichlet(0);
u.n[right] = dirichlet(0);
uf.n[bottom] =  dirichlet(0.);
uf.t[bottom] =  dirichlet(0.);

attribute {
	  double D;
}

double initial_mass;
event init (t = 0)
{
	 
	  foreach() {
		  if (y < tan(0.268)*(x-(1+0.338)))  {
			  f[] = 1.; //volfraction is 1 
			  f3[] = 0.;
	  	}
		  else if (x < ((y - tan(1.309)*(0.338))/(-tan(1.309))) && x > (( y/(-tan(1.309)))) && (y < 0.215) ) {
			  f[] = 1.;
			  f3[] = 1.;
		  }
			  else {
			  f[] = 0.;
		          f3[] = 0.;
		  }
		  foreach_dimension()
			  u.x[] = 0.;
		}
	
	/*
	fraction(f, y < 2 ||( y > 4 && y < 5 && x >3.5 && x < 4.5));
	foreach() {
		if(y > 4 && y < 5 && x > 3.5 && x < 4.5) {
			f3[] = f[];
		}
		foreach_dimension()
			u.x[]=0;
	}
*/
	//fraction(f3, x <= 0.338 && y < 0.215);
	boundary ({f,f3,u});
	initial_mass = statsf(f3).sum;

}


#include "view.h"
event logfile (i++){
	  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}


//Output .txt files for every timestep (x, y, f[])
static scalar * interfaces_save = NULL;

event vof (i++) {
	  boundary ({f, f3, uf});
	   
	    f3.inverse = false;
	      f3.gradient = minmod2;
	        f.tracers = {f3};

		  vof_advection({f}, i);

		    boundary({f, f3});

		      interfaces_save = interfaces;
		        interfaces = NULL;
}

event tracer_advection (i++) {
	  interfaces = interfaces_save;
}


#define D_L 0.0000001
mgstats mgd1;
event tracer_diffusion (i++) {
	  foreach()
		      f[] = clamp(f[], 0., 1.);
	    boundary ({f});

	      foreach()
		          f3[] = (f[] > F_ERR ? f3[]/f[] : 0.);
	        boundary({f3});
		  f3.D = D_L;
		    f3.inverse = false;
		      no_flux_diffusion (f3, f, dt);
		        foreach()
				    f3[] *= f[];
			  boundary({f3});
			

}




#define TRIANGLE (tan(0.268)*(x-(1+0.338+1.02)))
event makeslope (i++) {
		scalar tri[];
		fraction (tri, y <= TRIANGLE);
			foreach()
						foreach_dimension() 
									u.x[] -= u.x[]*tri[];


		face vector tri_face[];
 		 face_fraction (tri, tri_face);
 		 boundary ((scalar *){tri_face});
 		 foreach_face()
	   		 uf.x[] -= tri_face.x[]*uf.x[];
  		boundary ((scalar *){u, uf});  
			//	boundary((scalar *){u});

}






int time_step = 0.;
double shift_control = 1.;
event output_files  (t += 0.1; t < 2.2){ 
		
		
			foreach()
				f3[] = (f[] < 1e-6 ? 0 : f3[]);
			boundary({f3});
			
			double total_solute = statsf(f3).sum;
			//char s[80];
			//sprintf(s, "mass-%d.txt", t);
			//FILE * fp = fopen(s, "w");
			fprintf (fp0, "%.17g %.17g %.17g %.17g\n", t, initial_mass, total_solute,  total_solute - initial_mass);
			fflush (fp0);
	      
		        char one[80], two[80], three[80], four[80];
			scalar m[];		
			 foreach()
				     m[] = y < TRIANGLE;
			   boundary ({m});
			     	   
			//First two outputs show JPG snapshots for FOV 
			sprintf (one, "tracer-%d.png",time_step );
			output_ppm(f3, file = one, map = cool_warm, n  = 1028, box = {{0,0},{5,1}}, mask = m);
			sprintf (two, "density-%d.png",time_step );  
			output_ppm(rho, file = two, min = 1, max = 1400, map = cool_warm, linear = false, n = 640, box = {{0,0},{5,1}}, mask = m);

			view(psi = 0.268, tx = -0.55, ty = -0.05);
			box();
			char str[80];
			sprintf(str, "t = %g", t);
			draw_string(str, pos = 4, lw =3);
			draw_vof("f", filled = -1, fc = {1,1,1});
						squares("f3");
			save("movie.mp4");
			 time_step += 1.;
}



int add = 0;
event interface (t+= 0.1) {
	char fname[99],fname2[99];
	sprintf(fname, "water%d.csv", add);
	FILE * fp = fopen(fname, "w");
	output_facets(f,fp);
	fclose(fp);

}


int counter = 0;
#include "output_vtu_foreach.h"
event vtk_file (t += 0.1) {
		FILE * fp;
		char name[80];
		sprintf(name, "test_0.01_%d.vtu", counter);
		fp = fopen(name,"w");
		output_vtu_bin_foreach((scalar *) {rho,f}, (vector *) {u},N,fp,false);
		fclose(fp);
		counter += 1.;
}


/*

event adapt (i++)
{
	scalar omega[];
	vorticity(u,omega);
	 //adapt_wavelet((scalar *){cs,u,f},(double[]){ue,ue,ue,0.1},maxlevel,(maxlevel-4));
	adapt_wavelet((scalar *){cs,omega,f},(double[]){ue,ue,ue,0.01},maxlevel);
				        
}
*/
