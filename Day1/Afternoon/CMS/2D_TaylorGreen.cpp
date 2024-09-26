#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
using namespace std;
///----------------------------------------------------------------------------------
/// LB variables
const bool plot_vtk = true, central_moments = true;
const int np = 9;
const vector<int> cx = {0, 1, 0, -1,  0, 1, -1, -1,  1},
									cy = {0, 0, 1,  0, -1, 1,  1, -1, -1};
const vector<double> wf = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
const int nx = 10, ny = nx;
const double cs2 = 1./3., cs4 = 1./9., rho0 = 1., v0 = 0.01, Reynolds = 1000., kappa = 2.*M_PI/((double)nx), ni = v0*ny/Reynolds, T_ref = 1./(2.*kappa*kappa*ni), tau = ni*3.+0.5, omega = 1./tau, omega1 = 1.-omega;
const int nsteps = (int)(1*T_ref)+1, n_out = nsteps-1;
vector<double> f1(nx*ny*np, 0.), f2(nx*ny*np, 0.), u(nx*ny, 0.), v(nx*ny, 0.), rho(nx*ny, 0.), temp_pop(np, 0.);
double k1, k2, k3, k4, k5, k6, k7, k8;
double r0, r1, r2, r3, r4, r5, r6, r7, r8;
double U, V, R, ftemp, kinetic_energy, U2, V2, UV;
int newx, newy, id, idn, idnp;
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
	output_file << "X_COORDINATES " << nx << " float\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " float\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " float\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	/// Write density difference
	output_file << "SCALARS density float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
    {
      id = X*ny+Y;
      output_file << rho[id]-rho0<< "\n";
    }

	/// Write velocity
	output_file << "VECTORS velocity_vector float\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
    {
      id = X*ny+Y;
			output_file << u[id] << " " << v[id] << " 0\n";
    }

	/// Close file
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
			R = rho[id] = rho0/cs2*(1.-v0*v0/4./cs2*(cos(2.*kappa*x)+cos(2.*kappa*y)));
			U = u[id] = -v0*cos(kappa*x)*sin(kappa*y);
      V = v[id] = v0*sin(kappa*x)*cos(kappa*y);
			for(int k=0; k<np;k++)
        f1[id*np+k] = wf[k]*R*(1.+3.*(U*cx[k]+V*cy[k])+4.5*pow(U*cx[k]+V*cy[k],2)-1.5*(U*U+V*V));
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
void algoLB()
{
	int check = 0.;
	kinetic_energy = 0.;
  for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
      idnp = id*np;
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
			UV = U*V;

			r4 = f1[idnp+1]-f1[idnp+2]+f1[idnp+3]-f1[idnp+4];
			r5 = f1[idnp+5]-f1[idnp+6]+f1[idnp+7]-f1[idnp+8];

      if(central_moments==true)
      {
			k4 = r4-R*(U2-V2);
			k5 = r5-R*UV;
			k3 = 2*R*cs2;
	    k4 = omega1*k4;
			k5 = omega1*k5;
			k6 = 0.;
			k7 = 0.;
			k8 = R*cs4;
	    r0 = R;
			r1 = R*U;
			r2 = R*V;
			r3 = k3+R*(U2+V2);
			r4 = k4+R*(U2-V2);
			r5 = k5+R*UV;
			r6 = k6+2.*U*k5+0.5*V*(k3+k4)+R*U2*V;
			r7 = k7+2.*V*k5+0.5*U*(k3-k4)+R*U*V2;
			r8 = k8+2.*U*k7+2.*V*k6+4.*UV*k5-0.5*k4*(U2-V2)+0.5*k3*(U2+V2)+R*U2*V2;
		}
		else
		{
			r0 = R;
			r1 = R*U;
      r2 = R*V;
      r3 = R*(U2+V2+2.*cs2);
      r4 = r4*omega1+omega*R*(U2-V2);
      r5 = r5*omega1+omega*R*UV;
      r6 = R*V*(3*U2+cs2);
      r7 = R*U*(3*V2+cs2);
      r8 = R*(3*U2+1.)*(3.*V2+1.)*cs4;
			}
	    f1[idnp+0] = r0-r3+r8;
	    f1[idnp+1] = 0.5*(r1-r7-r8) + 0.25*(r3+r4);
	    f1[idnp+2] = 0.5*(r2-r6-r8) + 0.25*(r3-r4);
	    f1[idnp+3] = 0.25*(r3+r4) + 0.5*(-r1+r7-r8);
	    f1[idnp+4] = 0.25*(r3-r4) + 0.5*(-r2+r6-r8);
	    f1[idnp+5] = 0.25*(r5+r6+r7+r8);
	    f1[idnp+6] = 0.25*(r6-r5-r7+r8);
	    f1[idnp+7] = 0.25*(r5-r6-r7+r8);
	    f1[idnp+8] = 0.25*(r7-r6-r5+r8);

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
			kinetic_energy += R*(U2+V2);
      if(fabs(U)>1.)
        check = 1;
		}
	kinetic_energy /= (nx*ny*v0*v0);
  //return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_error()
{
	int time = nsteps-1;
	double u_an, v_an, vel_an, vel, num = 0., den = 0.;
  for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
    {
    	id = x*ny+y;
    	u_an = -v0*cos(kappa*x)*sin(kappa*y)*exp(-2.*kappa*kappa*time*ni);
    	v_an = v0*sin(kappa*x)*cos(kappa*y)*exp(-2.*kappa*kappa*time*ni);
    	vel_an = sqrt(u_an*u_an+v_an*v_an);
    	vel = sqrt(u[id]*u[id]+v[id]*v[id]);
    	num += pow(vel_an-vel,2);
    	den += pow(vel,2);
    }
   double error = sqrt(num/den);
   printf("Error = %e\n", error);

}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	//FILE *data = fopen("data.txt","wt");
	system("mkdir vtk_fluid");
	initial_state();
	clock_t c_start = std::clock();
	for(int i=0; i<nsteps; i++)
  {
		algoLB();
    f1 = f2;
		if(plot_vtk==true && i%n_out==0)
			write_fluid_vtk(i);
		//printf("t %lf of %lf\n", (double)i/T_ref, (double)nsteps/T_ref);
    //if(check_mach==1)
      //goto labelA;
  }
  compute_error();
  //labelA:
  //clock_t c_end = std::clock();
  //double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  //cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  //fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
