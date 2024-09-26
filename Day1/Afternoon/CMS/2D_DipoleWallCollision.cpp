#include <cstdio>
#include <cstring>
#include <cassert>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
///----------------------------------------------------------------------------------
using namespace std;
const bool plot_vtk = true;
const int nx = 512, ny = nx, np = 9;
const vector<int> cx = {0, 1, 0, -1, 0, 1, -1, -1, 1},
								  cy = {0, 0, 1, 0, -1, 1, 1, -1, -1},
								 opp = {0, 3, 4, 1 , 2, 7, 8, 5, 6};
const vector<double> wf = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
vector<double> f1(nx*ny*np, 0.), f2(nx*ny*np, 0.), u(nx*ny, 0.), v(nx*ny, 0.), rho(nx*ny, 0.), vorticity(nx*ny, 0.);
const double cs2 = 1./3., cs4 = cs2*cs2, rho0 = 1.;
double R, U, V, A, k0, k1, k2, k3, k4, k5, k6, k7, k8, ftemp, U2, V2, UV, r0, r1, r2, r3, r4, r5, r6, r7, r8, FX, FY, U2V2;
int check, newx, newy, id, idn;
double kinetic_energy, enstrophy, gradx_of_v, grady_of_u;
const double nx_phys = 2., we_phys = 299.56, T_ref_phys = 1, Scale_length = nx_phys/nx;
const double U_ref = 0.001, Scale_velocity = 1./U_ref, Scale_time = Scale_length/Scale_velocity, Scale_vorticity = 1./Scale_time, we = we_phys/Scale_vorticity;
const double Reynolds = 625., ni = nx/2.*U_ref/Reynolds, tau = ni*3+0.5, omega = 1./tau, omega1 = 1.-omega;
const double T_ref = T_ref_phys/Scale_time, R0 = 0.1/Scale_length;
const int nsteps = (int)(1*T_ref)+1, n_out = (int)(0.05*T_ref);
const double Scale_energy = pow(Scale_velocity,2)*pow(Scale_length,2), Scale_enstrophy = pow(Scale_velocity,2);
double Pxx, Pxy, Pyy, fneq, QPI, C;
vector<double> Feq(np,0.);
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	/// Create filename
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
	output_file << "X_COORDINATES " << nx << " double\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " double\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " double\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	output_file << "SCALARS vorticity float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			id = X*ny+Y;
			if(fabs(vorticity[id])<1e-16)
				vorticity[id] = 0.;
			output_file << vorticity[id]*Scale_vorticity << "\n";
		}

	/// Write velocity
	output_file << "VECTORS velocity_vector double\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			id = X*ny+Y;
			if(fabs(u[id])<1E-16)
				u[id] = 0.;
			if(fabs(v[id])<1E-16)
				v[id] = 0.;
			output_file << u[id] << " " << v[id] << " 0\n";
		}

	/// Close file
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	double i1 = nx/2, i2 = i1, j1 = ny/2+R0, j2 = ny/2-R0, R1, R2;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			rho[id] = rho0;
			R1 = sqrt((x-i1)*(x-i1)+(y-j1)*(y-j1));
			R2 = sqrt((x-i2)*(x-i2)+(y-j2)*(y-j2));
			U = u[id] = -0.5*we*(y-j1)*exp(-pow(R1/R0,2))+0.5*we*(y-j2)*exp(-pow(R2/R0,2));
      V = v[id] = 0.5*we*(x-i1)*exp(-pow(R1/R0,2))-0.5*we*(x-i2)*exp(-pow(R2/R0,2));
			for(int k=0; k<np;k++)
        f1[id*np+k] = f2[id*np+k] = wf[k]*rho[id]*(1.+3.*(U*cx[k]+V*cy[k])+4.5*pow(U*cx[k]+V*cy[k],2)-1.5*(U*U+V*V));
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
void walls()
{
	for(int x=0; x<nx; x++)
	{
		id = x*ny+0;
		for(int k=0; k<np; k++)
			f2[id*np+k] = f1[id*np+opp[k]];
		id = x*ny+ny-1;
		for(int k=0; k<np; k++)
			f2[id*np+k] = f1[id*np+opp[k]];
	}
	for(int y=0; y<ny; y++)
	{
		id = 0*ny+y;
		for(int k=0; k<np; k++)
			f2[id*np+k] = f1[id*np+opp[k]];
		id = (nx-1)*ny+y;
		for(int k=0; k<np; k++)
			f2[id*np+k] = f1[id*np+opp[k]];
	}
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_enstrophy()
{
	enstrophy = 0;
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			gradx_of_v = grady_of_u = 0;
			for(int k=0; k<np; k++)
			{
				newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
				idn = newx*ny+newy;
				gradx_of_v += 3.*wf[k]*v[idn]*cx[k];
				grady_of_u += 3.*wf[k]*u[idn]*cy[k];
			}
			vorticity[id] = gradx_of_v-grady_of_u;
			enstrophy += 0.5*pow(vorticity[id],2);
		}
/*	// Left
	int x = 0;
	for(int y=1; y<ny-1; y++)
	{
		id = x*ny+y;
		gradx_of_v = -1.5*v[id] + 2.*v[(x+1)*ny+y] - 0.5*v[(x+2)*ny+y];
		grady_of_u = -1.5*u[x*ny+(y-1)] + 1.5*u[x*ny+(y+1)];
	}
	enstrophy += 0.5*pow(gradx_of_v-grady_of_u,2);
	// Right
	x = nx-1;
	for(int y=1; y<ny-1; y++)
	{
		id = x*ny+y;
		gradx_of_v = 1.5*v[id] - 2.*v[(x-1)*ny+y] + 0.5*v[(x-2)*ny+y];
		grady_of_u = -1.5*u[x*ny+(y-1)] + 1.5*u[x*ny+(y+1)];
	}
	enstrophy += 0.5*pow(gradx_of_v-grady_of_u,2);
	//Bottom
	int y = 0;
	for(int x=1; x<nx-1; x++)
	{
		id = x*ny+y;
		gradx_of_v = -1.5*v[(x-1)*ny+y] + 1.5*v[(x+1)*ny+y];
		grady_of_u = -1.5*u[id] + 2.*u[x*ny+(y+1)] - 0.5*u[x*ny+(y+2)];
	}
	enstrophy += 0.5*pow(gradx_of_v-grady_of_u,2);
	//Top
	y = ny-1;
	for(int x=1; x<nx-1; x++)
	{
		id = x*ny+y;
		gradx_of_v = -1.5*v[(x-1)*ny+y] + 1.5*v[(x+1)*ny+y];
		grady_of_u = 1.5*u[id] - 2.*u[x*ny+(y-1)] + 0.5*u[x*ny+(y-2)];
	}
	enstrophy += 0.5*pow(gradx_of_v-grady_of_u,2);

	//BottomLeft
	x = 0;
	y = 0;
	id = x*ny+y;
	gradx_of_v = -1.5*v[id] + 2.*v[(x+1)*ny+y] - 0.5*v[(x+2)*ny+y];
	grady_of_u = -1.5*u[id] + 2.*u[x*ny+(y+1)] - 0.5*u[x*ny+(y+2)];
	enstrophy += 0.5*pow(gradx_of_v-grady_of_u,2);
	//BottomRight
	x = nx-1;
	y = 0;
	id = x*ny+y;
	gradx_of_v = 1.5*v[id] - 2.*v[(x-1)*ny+y] + 0.5*v[(x-2)*ny+y];
	grady_of_u = -1.5*u[id] + 2.*u[x*ny+(y+1)] - 0.5*u[x*ny+(y+2)];
	enstrophy += 0.5*pow(gradx_of_v-grady_of_u,2);
	//TopLeft
	x = 0;
	y = ny-1;
	id = x*ny+y;
	gradx_of_v = -1.5*v[id] + 2.*v[(x+1)*ny+y] - 0.5*v[(x+2)*ny+y];
	grady_of_u = 1.5*u[id] - 2.*u[x*ny+(y-1)] + 0.5*u[x*ny+(y-2)];
	enstrophy += 0.5*pow(gradx_of_v-grady_of_u,2);
	//TopRight
	x = nx-1;
	y = ny-1;
	id = x*ny+y;
	gradx_of_v = 1.5*v[id] - 2.*v[(x-1)*ny+y] + 0.5*v[(x-2)*ny+y];
	grady_of_u = 1.5*u[id] - 2.*u[x*ny+(y-1)] + 0.5*u[x*ny+(y-2)];
	enstrophy += 0.5*pow(gradx_of_v-grady_of_u,2);*/
}
///----------------------------------------------------------------------------------------------------------------------------------
int algoLB()
{
	int check = 0;
	kinetic_energy = 0.;
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			U = V = R = 0.;
			for(int k=0; k<np; k++)
			{
				ftemp = f1[id*np+k];
				R += ftemp;
				U += ftemp*cx[k];
				V += ftemp*cy[k];
			}
			rho[id] = R;
			U /= R;
			V /= R;
			u[id] = U;
			v[id] = V;
			U2 = U*U;
			V2 = V*V;
			U2V2 = U2*V2;
			UV = U*V;

			r4 = f1[id*np+1]-f1[id*np+2]+f1[id*np+3]-f1[id*np+4];
			r5 = f1[id*np+5]-f1[id*np+6]+f1[id*np+7]-f1[id*np+8];
			k4 = r4-R*(U2-V2);
			k5 = r5-R*UV;

	    k0 = R;
	    k1 = 0.;
	    k2 = 0.;
			k3 = 2.*R*cs2;
	    k4 = omega1*k4;
			k5 = omega1*k5;
			k6 = 0.;
			k7 = 0.;
			k8 = R*cs4;

	    r0 = k0;
			r1 = R*U;
			r2 = R*V;
			r3 = k3+R*(U2+V2);
			r4 = k4+R*(U2-V2);
			r5 = k5+R*UV;
			r6 = 2*U*k5+0.5*V*(k3+k4)+R*U2*V;
			r7 = 0.5*U*(k3-k4)+2*V*k5+R*U*V2;
			r8 = k8+0.5*k3*(U2+V2)-0.5*k4*(U2-V2)+R*U2V2+4*UV*k5;

	    f1[id*np+0] = r0-r3+r8;
	    f1[id*np+1] = 0.5*(r1-r7-r8) + 0.25*(r3+r4);
	    f1[id*np+2] = 0.5*(r2-r6-r8) + 0.25*(r3-r4);
	    f1[id*np+3] = 0.25*(r3+r4) + 0.5*(-r1+r7-r8);
	    f1[id*np+4] = 0.25*(r3-r4) + 0.5*(-r2+r6-r8);
	    f1[id*np+5] = 0.25*(r5+r6+r7+r8);
	    f1[id*np+6] = 0.25*(r6-r5-r7+r8);
	    f1[id*np+7] = 0.25*(r5-r6-r7+r8);
	    f1[id*np+8] = 0.25*(r7-r6-r5+r8);

			for(int k=0; k<np; k++)
			{
	    	newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
	      idn = newx*ny+newy;
				f2[idn*np+k] = f1[id*np+k];
			}
			kinetic_energy += 0.5*R*(U2+V2);
			if(fabs(U)>1)
				check = 1;
		}
	return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
void regularized()
{
	/// south & north
	for(int x=0; x<nx; x++)
	{
		int y = 0;
		id = x*ny+y;
		Pxx = Pxy = Pyy = 0.;
    U = u[id] = 0.;
  	V = v[id] = 0.;
  	R = rho[id] = rho0;
    C = -1.5*(U*U+V*V);
  	for(int k=0; k<np; k++)
  	{
    	A = U*cx[k]+V*cy[k];
 	     Feq[k] = wf[k]*R*(1.+3.*A+4.5*A*A+C);
  		fneq = f1[id*np+k]-Feq[k];
    	Pxx += fneq*cx[k]*cx[k];
    	Pxy += fneq*cx[k]*cy[k];
  		Pyy += fneq*cy[k]*cy[k];
    }
    for(int k=0; k<np; k++)
  	{
  		QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
    	f2[id*np+k] = Feq[k]+4.5*wf[k]*QPI;
    }
		y = ny-1;
		id = x*ny+y;
		Pxx = Pxy = Pyy = 0.;
  	U = u[id] = 0.;
  	V = v[id] = 0.;
    R = rho[id] = rho0;
  	C = -1.5*(U*U+V*V);
  	for(int k=0; k<np; k++)
    {
  		A = U*cx[k]+V*cy[k];
 	    Feq[k] = wf[k]*R*(1.+3.*A+4.5*A*A+C);
    	fneq = f1[id*np+k]-Feq[k];
    	Pxx += fneq*cx[k]*cx[k];
    	Pxy += fneq*cx[k]*cy[k];
  		Pyy += fneq*cy[k]*cy[k];
  	}
    for(int k=0; k<np; k++)
  	{
  		QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
    	f2[id*np+k] = Feq[k]+4.5*wf[k]*QPI;
    }
	}

	/// west & east
	for(int y=0; y<nx; y++)
	{
		int x = 0;
		id = x*ny+y;
		Pxx = Pxy = Pyy = 0.;
    U = u[id] = 0.;
  	V = v[id] = 0.;
  	R = rho[id] = rho0;
    C = -1.5*(U*U+V*V);
  	for(int k=0; k<np; k++)
  	{
    	A = U*cx[k]+V*cy[k];
 	    Feq[k] = wf[k]*R*(1.+3.*A+4.5*A*A+C);
  		fneq = f1[id*np+k]-Feq[k];
    	Pxx += fneq*cx[k]*cx[k];
    	Pxy += fneq*cx[k]*cy[k];
  		Pyy += fneq*cy[k]*cy[k];
    }
    for(int k=0; k<np; k++)
  	{
  		QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
    	f2[id*np+k] = Feq[k]+4.5*wf[k]*QPI;
    }
		x = nx-1;
		id = x*ny+y;
		Pxx = Pxy = Pyy = 0.;
  	U = u[id] = 0.;
  	V = v[id] = 0.;
    R = rho[id] = rho0;
  	C = -1.5*(U*U+V*V);
  	for(int k=0; k<np; k++)
    {
  		A = U*cx[k]+V*cy[k];
 	    Feq[k] = wf[k]*R*(1.+3.*A+4.5*A*A+C);
    	fneq = f1[id*np+k]-Feq[k];
    	Pxx += fneq*cx[k]*cx[k];
    	Pxy += fneq*cx[k]*cy[k];
  		Pyy += fneq*cy[k]*cy[k];
  	}
    for(int k=0; k<np; k++)
  	{
  		QPI = Pxx*(cx[k]*cx[k]-cs2)+2*Pxy*cx[k]*cy[k]+Pyy*(cy[k]*cy[k]-cs2);
    	f2[id*np+k] = Feq[k]+4.5*wf[k]*QPI;
    }
	}
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  system("mkdir vtk_fluid");
	FILE *data_output = fopen("data.txt","wt");
	int check_mach = 0;
	initial_state();
	//printf("%lf %lf %lf %lf %lf %lf\n", Scale_length, Scale_time, Scale_energy, Scale_enstrophy, we, T_ref);
  for(int i=0; i<nsteps; i++)
  {
    check_mach = algoLB();
		//walls();
		regularized();
		compute_enstrophy();
		f1 = f2;
		fprintf(data_output,"%lf    %e   %e\n", i*Scale_time, kinetic_energy*Scale_energy, enstrophy*Scale_enstrophy);
		if(i%n_out==0)
			printf("step = %lf of %lf. Energy = %lf, Enstrophy = %lf\n", i*Scale_time, nsteps*Scale_time, kinetic_energy*Scale_energy, enstrophy*Scale_enstrophy);
		if(i%n_out==0 && plot_vtk==true)
			write_fluid_vtk(i);
		if(check_mach==1)
      goto labelA;
  }
  labelA:
  fclose(data_output);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
