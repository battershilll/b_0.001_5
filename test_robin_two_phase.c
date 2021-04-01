#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
//#include "tracer.h"

//# define rho(f) (clamp(f,0.,1.)*(clamp(s[],0.,1.)*(rho3 - rho1) + \
//		rho1 - rho2) + rho2)
scalar f3[];
double rho3 = 1400; //density of dense salt water
double mu3 = 0.1; //density of salt water at 20degrees

//scalar * tracers = {f3};

//We “overload” the default by defining the mu() macro before including the code for two phases.
/*
#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(clamp(f3[],0.,1.)*(rho3 - rho1) + rho1 -rho2) + rho2) //void fraction (density) is a function of the volume fraction. can add on a third phase here
#endif
#ifndef mu
# define mu(f)  (clamp(f,0.,1.)*(clamp(f3[],0.,1.)*(mu3 - mu1) + mu1 - mu2) + mu2) //viscoisty is a function of the volumne fraction
#endif
*/

#define MAXLEVEL 11

#define robin(a,b,c) ((dirichlet ((c)*Delta/(2*(b) + (a)*Delta))) + ((neumann (0))* ((2*(b) - (a)*Delta)/(2*(b) + (a)*Delta) + 1.)))
//#include "conserving.h"

/*
#ifndef diff
# define diff(f3) (clamp(f3,0.,1)*0.1)
#endif*/


//#include "advection_Q.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h" //updated conserving
//#include "navier-stokes/conserving.h"
#include "mixtures.h"
//#include "tension.h"
//#include "elementary_body.h"

#include "utils.h"

//scalar * tracers = {f3}; //f3 is a tracer for the denstiy of the third phase


//rho3 = 1500;
double ue = 1e-3;
FILE * fp;

#define io (0.338)
#define MUR 54
/**
  The domain is eight units long, centered vertically. */

int main() {
	  L0 = 8.;
	    origin (0, 0);
	    rho1 = 998; //water
	    rho2 = 1.27; //air
	    //double rho3 = 1400;
	    mu1 = 1.002e-3; //water
            mu2 = 1.002e-3/MUR; //air	    
	    //mu3 = 0.0001:;
	    //init_grid(2048);
	    //f.sigma = 73e-3;
	    //G.x = 0.;
	    //G.y = -9.81;
	    //init_grid (16.);
	    init_grid (1 << MAXLEVEL);
	    fp=fopen("tracers_cons_final.dat", "w");
 	    const face vector av[] = {9.81*sin(0.268), -9.81*cos(0.268), 0};
	    a = av;
	    //DT = 0.001;
	   // TOLERANCE = 1e-;
	    //CFL = 0.4;
	    run(); 
}



//uf.n[bottom] =  dirichlet(0.);
//uf.t[bottom] =  dirichlet(0.);


//DEFAULT BOUNDARY CONDIITONS ARE NEUMANN (free slip)
//u.n[bottom] =  robin(0.5,0.5,0.); //robin BC
u.t[bottom] = robin(1.,0.0,0.); //robin BC - 0.015 should be the height of the flow. As B --> 0.015, the condition tends towards no slip

u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);

//f3[bottom] = neumann(0.); //free slip for the tracer
//u.t[bottom] =  dirichlet(0. * (t < 0.7));



attribute {
	  double D;
}


event init (t = 0)
{
	 
	  //refine(y < 2. && level < MAXLEVEL);
	  foreach() {
		  if (y < tan(0.268)*(x-(1+io))) { //with water
			  f[] = 1.;
			  //f3[] = 0.;
		  }
		  else if (x < ((y - tan(1.309)*(io))/(-tan(1.309))) && y < 0.185 ) { //accounts for the extra volume
			  f[] = 1.;
			  //f3[] = 1.;
		  }
			  else {
			  f[] = 0.;
		          //f3[] = 0.;
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
	boundary ({f,u});

}


#include "view.h"
event logfile (i++){
	  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}


//vof event is in conserving.h
/*
#define D_L 0.0001
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
*/

// Bottom boundary 
scalar tri[];
#define TRIANGLE (tan(0.268)*(x-(1+io+1.02)))
event makeslope (i++) {
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




double pe1a = 0.00;
event energy_output_full (t+=0.01) {
	static FILE * fpET = fopen("energy_total.txt", "w");
	double ke1 = 0.000000000;
	double pe1 = 0.000000000;
	//double pe1a = 0.00000000;

	if (t==0.0) {

		foreach(reduction(+:pe1a)) {
			pe1a += rho(f[])*dv()*9.81*(-sin(0.268)*x+cos(0.268)*y);
		}
	}

	foreach(reduction(+:ke1) reduction(+:pe1)) {
 		double r = rho(f[]);
 		ke1 += r*sq(norm(u))/2.*dv();
 		pe1 += (r*dv()*9.81*(-sin(0.268)*x+cos(0.268)*y));
		//pe1 += r*dv()*9.81*cos(0.268)*y+r*dv()*9.81*sin(0.268)*x;
	}
	fprintf (fpET, "%g %g %g\n", t, ke1, pe1 -  pe1a);
	fflush(fpET);
	

	static FILE * fpEI = fopen("energy_init.txt", "w");
	double ke2 = 0.000000000;
	double pe2 = 0.000000000;

	foreach(reduction(+:ke2) reduction(+:pe2)) {
		if (x <= 1.338) {
 			double r = rho(f[]);
 			ke2 += r*sq(norm(u))/2.*dv();
 			pe2 += (r*dv()*9.81*(-sin(0.268)*x+cos(0.268)*y));
		}
		//pe1 += r*dv()*9.81*cos(0.268)*y+r*dv()*9.81*sin(0.268)*x;
	}
	fprintf (fpEI, "%g %g %g\n", t, ke2, pe2);
	fflush(fpEI);


	static FILE * fpEE = fopen("energy_end.txt", "w");
	double ke3 = 0.000000000;
	double pe3 = 0.000000000;

	foreach(reduction(+:ke3) reduction(+:pe3)) {
		if (x > 7.00) {
 		double r = rho(f[]);
 		ke3 += r*sq(norm(u))/2.*dv();
 		pe3 += (r*dv()*9.81*(-sin(0.268)*x+cos(0.268)*y));
		}
		//pe1 += r*dv()*9.81*cos(0.268)*y+r*dv()*9.81*sin(0.268)*x;
	}
	fprintf (fpEE, "%g %g %g\n", t, ke3, pe3);
	fflush(fpEE);


	  


}



/*
event graphs (i++) {
	static FILE * fptot = fopen("budget_total.dat", "w");
    	//static FILE * fpair = fopen("budgetAir.dat", "w");
	double ke = 0., gpe = 0.;
	//double keAir = 0., gpeAir = 0.;
	foreach(reduction(+:ke) reduction(+:gpe))  {
		double norm2 = 0.;
		foreach_dimension()
			norm2 += sq(u.x[]);
		ke += rho(f[])*norm2*dv();
		//keAir += rho[]*norm2*(1.0-f[])*dv();
		gpe += rho(f[])*9.81*cos(0.268)*y*dv() + rho(f[])*9.81*sin(0.268)*x*dv();
		//gpeAir += rho2*g_*y*(1.0-f[])*dv();
	}
	double rates[2];
	dissipation_rate(rates);
	double dissWater = rates[0];
	double dissAir   = rates[1];
	if (i == 0) {
	fprintf (fpwater, "t ke gpe dissipation\n");
	//fprintf (fpair, "t ke gpe dissipation\n");
	}
	fprintf (fpwater, "%g %g %g %g\n",
			t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, dissWater);
	//fprintf (fpair, "%g %g %g %g\n",
	//		t/(k_/sqrt(g_*k_)), keAir/2., gpeAir + 0.125, dissAir);
	fprintf (ferr, "%g %g %g %g\n",
			t/(k_/sqrt(g_*k_)), ke/2., gpe + 0.125, dissWater);
}
*/






int time_step = 0.;
double shift_control = 1.;
event output_files  (t += 0.01; t < 4.5) { 
		
		
			
			char one[80], two[80], three[80], four[80];
				
   			  foreach() {
				      if (tri[] > 0.5)
					            tri[] = -1;
				        }
		     	   
			//First two outputs show JPG snapshots for FOV 
			sprintf (one, "tracer-%d.png",time_step );
			output_ppm(f, file = one, map = cool_warm, n = 2048, linear = true, box = {{0,0},{4.5,1}}, mask = tri);
			sprintf (two, "density-%d.png",time_step );  
			output_ppm(rho, file = two, min = 1, max = 1400, n = 2048, linear = true, box = {{0,0},{4.5,1}}, mask = tri);
			//sprintf (three, "velocity_profile-%d.txt", time_step);
			//output_field(u.x, file = three, linear = true, box = {{1.238,0},{1.288 +0.01,0.05}});
	
			view(fov = 20, psi = 0.268, tx = -0.55, ty = -0.05);
			box();
			char str[80];
			sprintf(str, "t = %g", t);
			draw_string(str, pos = 4, lw =3);
			draw_vof("f", filled = -1, fc = {1,1,1});
						squares("f");
			save("movie.mp4");
			
			//We also output a vertical cross section of the velocity compoentns 

			//FILE*fp1 = fopen(three, "w");
			sprintf (three, "vprof-%d",time_step);
			FILE * fp1 = fopen (three, "w");
			//double step = (L0)/(1 << (lmax));
			for (double y = 0.; y <= 0.05; y += 8/(pow(2,MAXLEVEL))){
				fprintf (fp1, "%g %g %g %g\n", y,interpolate (u.x, 1.238, y),interpolate (u.y, 1.238, y),interpolate(f, 1.238,y));} //Outputs y, u.x, u.y, f
			fclose (fp1);

			time_step += 1.;
}

//Output interface

int add = 0;
event interface (t+= 0.01) {
		char fname[99];
		sprintf(fname, "water_elevation%d.txt", add);
		FILE * fp = fopen(fname, "w");
		output_facets(f,fp);
		fclose(fp);
		add += 1.;

}

//Output vtu files

int counter = 0;
#include "output_vtu_foreach.h"
event vtk_file (t += 0.05) {
			FILE * fp;
			char name[80];
			sprintf(name, "test_0.01_%d.vtu", counter);
			fp = fopen(name,"w");
			output_vtu_bin_foreach((scalar *) {rho,f}, (vector *) {u},N,fp,false);
			fclose(fp);
			counter += 1.;
}






event adapt (i++)
{
	scalar omega[];
	vorticity(u,omega);
	 //adapt_wavelet((scalar *){cs,u,f},(double[]){ue,ue,ue,0.1},maxlevel,(maxlevel-4));
	adapt_wavelet((scalar *){omega,f},(double[]){ue,ue,0.01},MAXLEVEL);
				        
}

