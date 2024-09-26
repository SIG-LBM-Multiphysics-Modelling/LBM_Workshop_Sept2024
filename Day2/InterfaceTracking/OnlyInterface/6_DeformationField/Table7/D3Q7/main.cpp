#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
const bool plot_vtk = true;
const int nx = 100, ny = nx, nz = nx, np = 7;
//          						0  1   2  3   4  5   6  7   8   9  10 11  12  13  14 15  16  17  18
vector<const int> cx = {0, 1, -1, 0,  0, 0,  0},
									cy = {0, 0,  0, 1, -1, 0,  0},
									cz = {0, 0,  0, 0,  0, 1, -1};
const double cs2 = 1./4., cs4 = cs2*cs2, d = (double)(nx), U_ref = 0.02, T_ref = d/U_ref, Pe = 200, Ch = 0.03, xi = 3;
vector<const double> wf = {1./4., 1/8., 1/8., 1/8., 1/8., 1/8., 1/8.};
const int nsteps = (int)(1.*T_ref+1), n_out = (int)(T_ref/4);
const double rhoL = 1., rhoH = rhoL, sigma = 1E-5;
vector<double> f1(nx*ny*nz*np,0.), f2(nx*ny*nz*np,0.), rho(nx*ny*nz,0.), press(nx*ny*nz,0.), u(nx*ny*nz,0.), v(nx*ny*nz,0.), w(nx*ny*nz,0.), temp_pop(np,0.);
double feq, tau, omega, omega1, ni, Fpx, Fpy, Fpz, Fsx, Fsy, Fsz, Fmx, Fmy, Fmz, mu, Press, Press0, U0, V0, W0;
double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18;
double r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18;
double R, U, V, W, ftemp, CX, CY, CZ, A, B, C, U2, V2, W2, Fx, Fy, Fz, U3, V3, W3, U4, V4, W4, P1;
int newx, newy, newz, id, idn;
// PHASE
double k1_g, k2_g, k3_g, k4_g, k5_g, k6_g, k7_g, k8_g, k9_g, k10_g, k11_g, k12_g, k13_g, k14_g, k15_g, k16_g, k17_g, k18_g;
const double PhiH = 1., PhiL = 0., Phi0 = 0.5*(PhiH+PhiL), beta = 12.*sigma/xi, kappa = 3.*sigma*xi/2.;
const double M = 0.001, tau_phase = M/cs2+0.5, omega_phase = 1./tau_phase, omega_phase1 = 1.-omega_phase;
double grad_phix, grad_phiy, grad_phiz, laplPhi, Phi, gtemp, P, Fx_phase, Fy_phase, Fz_phase, Nx, Ny, Nz, Phii, thisPhi, Phi_prev;
vector<double> g1(nx*ny*nz*np,0.), g2(nx*ny*nz*np,0.), phase(nx*ny*nz,0.), phase_old(nx*ny*nz,0.), temp_pop_phase(np,0.), phase0(nx*ny*nz,0.);
double error;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;
	output_file.open(output_filename.str().c_str());

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " " << nz << " " << "\n";
	output_file << "X_COORDINATES " << nx << " float\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " float\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << nz  << " float\n";
	for(int z = 0; z < nz ; ++z)
		output_file << z  << " ";
	output_file << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) * (nz)  << "\n";

	output_file << "SCALARS phase float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Z = 0; Z < nz; ++Z)
		for(int Y = 0; Y < ny ; ++Y)
			for(int X = 0; X < nx; ++X)
			{
				id = (X*ny+Y)*nz+Z;
				if(phase[id]>1E-12)
					output_file << phase[id] << "\n";
				else
					output_file << "0 \n";
			}
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	double X, Y, Z, dist;
	int radius = nx/5;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
    	for(int z=0; z<nz; z++)
			{
				X = double(x) / double(nx-1);
				Y = double(y) / double(ny-1);
				Z = double(z) / double(nz-1);
				id = (x*ny+y)*nz+z;
				dist = sqrt((x-nx/2)*(x-nx/2)+(y-ny/2)*(y-ny/2)+(z-nz/2)*(z-nz/2));
				phase[id] = 0.5*(PhiH+PhiL)+0.5*(PhiH-PhiL)*tanh(2./xi*(radius-dist));
				// if( (x-3*nx/10)*(x-3*nx/10)+(y-3*ny/10)*(y-3*ny/10)+(z-nz/2)*(z-nz/2)<radius*radius )
				// 	phase[id] = PhiH;
			  // else
				// 	phase[id] = PhiL;
				U = 0.5*U_ref*(sin(4*M_PI*(X-0.5)) * sin(4*M_PI*(Y-0.5)) + cos(4*M_PI*(Z-0.5)) * cos(4*M_PI*(X-0.5)));
      	V = 0.5*U_ref*(sin(4*M_PI*(Y-0.5)) * sin(4*M_PI*(Z-0.5)) + cos(4*M_PI*(X-0.5)) * cos(4*M_PI*(Y-0.5)));
	      W = 0.5*U_ref*(sin(4*M_PI*(Z-0.5)) * sin(4*M_PI*(X-0.5)) + cos(4*M_PI*(Y-0.5)) * cos(4*M_PI*(Z-0.5)));
				C = -1.5*(U*U+V*V+W*W);
				for(int k=0; k<np; k++)
				{
	        A = U*cx[k]+V*cy[k]+W*cz[k];
	        B = 4.5*A*A;
					g1[id*np+k]= wf[k]*phase[id]*(1.+A/cs2);
				}
			}
	phase0 = phase;
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_grad_phi(int x, int y, int z)
{
  grad_phix = grad_phiy = grad_phiz = laplPhi = 0.;
  thisPhi = phase_old[id];
	for(int k=1; k<np; ++k)
	{
		newx = x+cx[k];
		newy = y+cy[k];
		newz = z+cz[k];
		if(x==0 || x==nx-1)
			newx = (newx+nx)%nx;
		if(y==0 || y==ny-1)
			newy = (newy+ny)%ny;
		if(z==0 || z==nz-1)
			newz = (newz+nz)%nz;
		idn = (newx*ny+newy)*nz+newz;
		Phii = phase_old[idn];
		grad_phix += Phii*wf[k]*cx[k];
		grad_phiy += Phii*wf[k]*cy[k];
		grad_phiz += Phii*wf[k]*cz[k];
    laplPhi += wf[k]*(Phii-thisPhi);
	}
	grad_phix /= cs2;
	grad_phiy /= cs2;
	grad_phiz /= cs2;
  laplPhi *= 2./cs2;
}
///----------------------------------------------------------------------------------------------------------------------------------
int algo_LB(int time)
{
	int hh = 0;
	double X, Y, Z;
	phase_old = phase;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
    	for(int z=0; z<nz; z++)
			{
				id = (x*ny+y)*nz+z;
				//compute_grad_phi(x, y, z);
				Phi = 0.;
				for(int k=0; k<np; ++k)
				{
					temp_pop_phase[k] = gtemp = g1[id*np+k];
					Phi += gtemp;
				}
				phase[id] = Phi;
				//rho[id] = R = rhoL+(Phi_prev-PhiL)/(PhiH-PhiL)*(rhoH-rhoL);
				X = double(x) / double(nx-1);
				Y = double(y) / double(ny-1);
				Z = double(z) / double(nz-1);
				U = 0.5*U_ref*(sin(4*M_PI*(X-0.5)) * sin(4*M_PI*(Y-0.5)) + cos(4*M_PI*(Z-0.5)) * cos(4*M_PI*(X-0.5))) * cos(M_PI*time/T_ref);
      	V = 0.5*U_ref*(sin(4*M_PI*(Y-0.5)) * sin(4*M_PI*(Z-0.5)) + cos(4*M_PI*(X-0.5)) * cos(4*M_PI*(Y-0.5))) * cos(M_PI*time/T_ref);
	      W = 0.5*U_ref*(sin(4*M_PI*(Z-0.5)) * sin(4*M_PI*(X-0.5)) + cos(4*M_PI*(Y-0.5)) * cos(4*M_PI*(Z-0.5))) * cos(M_PI*time/T_ref);
				U2 = U*U;
				V2 = V*V;
				W2 = W*W;
				if(fabs(U)>1.)
          hh = 1;

				r1 = temp_pop_phase[1]-temp_pop_phase[2];
				r2 = temp_pop_phase[3]-temp_pop_phase[4];
				r3 = temp_pop_phase[5]-temp_pop_phase[6];
				k1_g = r1-Phi*U;
				k2_g = r2-Phi*V;
				k3_g = r3-Phi*W;
				Nx = -k1_g/(sqrt(pow(k1_g,2)+pow(k2_g,2)+pow(k3_g,2))+1e-12);
				Ny = -k2_g/(sqrt(pow(k1_g,2)+pow(k2_g,2)+pow(k3_g,2))+1e-12);
				Nz = -k3_g/(sqrt(pow(k1_g,2)+pow(k2_g,2)+pow(k3_g,2))+1e-12);


				for(int k=0; k<np; ++k)
				{
					A = U*cx[k]+V*cy[k]+W*cz[k];
					g1[id*np+k] = omega_phase1*g1[id*np+k] + omega_phase* ( wf[k]*Phi*(1.+A/cs2) + wf[k]*M/cs2*4*Phi*(1-Phi)/xi*(cx[k]*Nx+cy[k]*Ny+cz[k]*Nz) );
					newx = x+cx[k];
					newy = y+cy[k];
					newz = z+cz[k];
					if(x==0 || x==nx-1)
						newx = (newx+nx)%nx;
					if(y==0 || y==ny-1)
						newy = (newy+ny)%ny;
					if(z==0 || z==nz-1)
						newz = (newz+nz)%nz;
					idn = (newx*ny+newy)*nz+newz;
					g2[idn*np+k] = g1[id*np+k];
				}
			}
  return hh;
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_error()
{
	double num = 0., den = 0.;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
    	for(int z=0; z<nz; z++)
			{
				id = (x*ny+y)*nz+z;
				num += pow(phase0[id]-phase[id],2);
				den += pow(phase0[id],2);
			}
	error = sqrt(num/den);
	printf("Error = %lf. \n", error);
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  system("mkdir vtk_fluid");
  FILE *data = fopen("data.txt","wt");
	initial_state();
	int check_mach = 0;
	clock_t c_start = clock();
	printf("%lf %lf\n", rhoH, rhoL);
	for(int i=0; i<nsteps; i++)
  {
    check_mach = algo_LB(i);
		g1 = g2;
    //fprintf(data,"%lf    %e    %e\n", (double)i/T_ref, kinetic_energy, enstrophy);
    if(check_mach==1)
      goto labelA;
		if(plot_vtk==true && i%n_out==0)
			write_fluid_vtk(i);
   	if(i%1==0)
      printf("Iteration %lf of %lf. \n", (double)i/T_ref, (double)nsteps/T_ref);
  }
compute_error();
	labelA:
  clock_t c_end = clock();
	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";
	compute_error();
	fprintf(data,"Error=%lf\n", error);
	fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
