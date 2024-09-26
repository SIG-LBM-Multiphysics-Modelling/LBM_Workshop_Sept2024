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
///----------------------------------------------------------------------------------
const bool plot_vtk = true;
const int nx = 64, ny = nx, nz = nx, np = 19;
vector<const int> cx = {0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1, 0,  0,  0,  0},
									cy = {0, 0,  0, 1, -1, 0,  0, 1, -1, -1,  1, 0,  0,  0,  0, 1, -1,  1, -1},
									cz = {0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1, -1, -1,  1, 1, -1, -1,  1};
const double cs = 1./sqrt(3.), cs2 = 1./3., cs4 = cs2*cs2, rho0 = 1.;
const double Ma = 0.1, v0 = Ma*cs, Reynolds = 1600., ni = v0*(nx-1)/Reynolds, tau = ni/cs2+0.5, omega = 1./tau, omega1 = 1.-omega;
const double T_ref = ((double)nx-1)/(v0);
const int nsteps = (int)(2*T_ref)+1, n_out = (int)(T_ref/1)+1;
vector<const double> wf = {1/3., 1/18., 1/18., 1/18., 1/18., 1/18., 1/18., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36.};
vector<double> f1(nx*ny*nz*np,0.), f2(nx*ny*nz*np,0.), u(nx*ny*nz,0.), v(nx*ny*nz,0.), w(nx*ny*nz,0.);
double R, U, V, W, ftemp, feq, kinetic_energy, kinetic_energy_0, U2, V2, W2;
float k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18;
float r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18;
int newx, newy, newz, check, id, idn;
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

	/// Write velocity
	output_file << "VECTORS velocity_vector float\n";
	for(int x=0; x<nx; ++x)
		for(int y=0; y<ny; ++y)
			for(int z=0; z<nz; ++z)
			{
				id = (x*ny+y)*nz+z;
				output_file << u[id] << " " << v[id] << " " << w[id] << "\n";
			}
	/// Close file
	output_file.close();
}
///---------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	kinetic_energy_0 = 0.;
	double II, JJ, ZZ;
	for(int x=0; x<nx; ++x)
		for(int y=0; y<ny; ++y)
			for(int z=0; z<nz; ++z)
			{
				id = (x*ny+y)*nz+z;
				II = 2*M_PI*(double)x/(double)(nx-1);
      	JJ = 2*M_PI*(double)y/(double)(ny-1);
      	ZZ = 2*M_PI*(double)z/(double)(nz-1);
				R = rho0;
				U = u[id] = v0*cos(II)*sin(JJ)*sin(ZZ);
        V = v[id] = -v0*sin(II)*cos(JJ)*sin(ZZ)*0.5;
        W = w[id] = -v0*sin(II)*sin(JJ)*cos(ZZ)*0.5;
				for(int k=0; k<np;k++)
					f1[id*np+k] = wf[k] * R * (1. + 1./cs2*(U*cx[k]+V*cy[k]+W*cz[k]) + 0.5/cs2/cs2*pow(U*cx[k]+V*cy[k]+W*cz[k],2) - 0.5/cs2*(U*U+V*V+W*W));
				kinetic_energy_0 += R*(U*U+V*V+W*W);
			}
}
///----------------------------------------------------------------------------------------------------------------------------------
int algo_LB()
{
	check = 0.;
	kinetic_energy = 0.;
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
			for(int z=0; z<nz; z++)
			{
				id = (x*ny+y)*nz+z;
				U = V = W = R = 0.;
				for(int k=0; k<np; k++)
				{
					ftemp = f1[id*np+k];
					R += ftemp;
					U += ftemp*cx[k];
					V += ftemp*cy[k];
					W += ftemp*cz[k];
				}
				U /= R;
				V /= R;
				W /= R;
				u[id] = U;
				v[id] = V;
				w[id] = W;

				r5 = f1[id*np+1]+f1[id*np+2]-f1[id*np+3]-f1[id*np+4]+f1[id*np+11]+f1[id*np+12]+f1[id*np+13]+f1[id*np+14]-f1[id*np+15]-f1[id*np+16]-f1[id*np+17]-f1[id*np+18];
				r6 = f1[id*np+3]+f1[id*np+4]-f1[id*np+5]-f1[id*np+6]+f1[id*np+7]+f1[id*np+8]+f1[id*np+9]+f1[id*np+10]-f1[id*np+11]-f1[id*np+12]-f1[id*np+13]-f1[id*np+14];
				r7 = f1[id*np+7]+f1[id*np+8]-f1[id*np+9]-f1[id*np+10];
				r8 = f1[id*np+11]+f1[id*np+12]-f1[id*np+13]-f1[id*np+14];
				r9 = f1[id*np+15]+f1[id*np+16]-f1[id*np+17]-f1[id*np+18];
				k5 = r5-R*(U2-V2);
				k6 = r6-R*(V2-W2);
				k7 = r7-R*U*V;
				k8 = r8-R*U*W;
				k9 = r9-R*V*W;
				///collide moments
				k0 = R;
				k1 = 0.;
				k2 = 0.;
				k3 = 0.;
				k4 = R;
				k5 = omega1*k5;
				k6 = omega1*k6;
				k7 = omega1*k7;
				k8 = omega1*k8;
				k9 = omega1*k9;
				k10 = 0.;
				k11 = 0.;
				k12 = 0.;
				k13 = 0.;
				k14 = 0.;
				k15 = 0.;
				k16 = R*cs4;
				k17 = R*cs4;
				k18 = R*cs4;

				r0 = R;
				r1 = k1+R*U;
				r2 = k2+R*V;
				r3 = k3+R*W;
				r4 = k4+2.*U*k1+2.*V*k2+2.*W*k3+R*(U2+V2+W2);
				r5 = k5+2.*U*k1-2.*V*k2+R*(U2-V2);
				r6 = k6+2.*V*k2-2.*W*k3+R*(V2-W2);
				r7 = k7+U*k2+V*k1+R*U*V;
				r8 = k8+U*k3+W*k1+R*U*W;
				r9 = k9+V*k3+W*k2+R*V*W;
				r10 = k10+2.*U*k7+U2*k2+2.*U*V*k1+V*cs2*k4+2.*V*cs2*k5+V*cs2*k6+R*U2*V;
				r11 = k11+2.*V*k7+V2*k1+2.*U*V*k2+U*cs2*k4-U*cs2*k5+U*cs2*k6+R*U*V2;
				r12 = k12+2.*U*k8+U2*k3+2.*U*W*k1+W*cs2*k4+2.*W*cs2*k5+W*cs2*k6+R*U2*W;
				r13 = k13+2.*W*k8+W2*k1+2.*U*W*k3+U*cs2*k4-U*cs2*k5-2.*U*cs2*k6+R*U*W2;
				r14 = k14+2.*V*k9+V2*k3+2.*V*W*k2+W*cs2*k4-W*cs2*k5+W*cs2*k6+R*V2*W;
				r15 = k15+2.*W*k9+W2*k2+2.*V*W*k3+V*cs2*k4-V*cs2*k5-2.*V*cs2*k6+R*V*W2;
				r16 = k16+k4*(cs2*U2+cs2*V2)+k6*(cs2*U2+cs2*V2)-k5*(cs2*U2-2.*cs2*V2)+2.*U*k11+2.*V*k10+R*U2*V2+4*U*V*k7+2.*U*V2*k1+2.*U2*V*k2;
				r17 = k17+k4*(cs2*U2+cs2*W2)-k5*(cs2*U2-2.*cs2*W2)-k6*(2.*cs2*U2-cs2*W2)+2.*U*k13+2.*W*k12+R*U2*W2+4*U*W*k8+2.*U*W2*k1+2.*U2*W*k3;
				r18 = k18+k4*(cs2*V2+cs2*W2)-k5*(cs2*V2+cs2*W2)-k6*(2.*cs2*V2-cs2*W2)+2.*V*k15+2.*W*k14+R*V2*W2+4*V*W*k9+2.*V*W2*k2+2.*V2*W*k3;

				f1[id*np+0] = r0-r4+r16+r17+r18;
				f1[id*np+1] = (r4+2.*r5+r6)*0.5*cs2-(-r1+r11+r13+r16+r17)*0.5;
				f1[id*np+2] = (r4+2.*r5+r6)*0.5*cs2+(r11+r13-r1-r16-r17)*0.5;
				f1[id*np+3] = (r4-r5+r6)*0.5*cs2-(-r2+r10+r15+r16+r18)*0.5;
				f1[id*np+4] = (r4-r5+r6)*0.5*cs2-(r2-r10-r15+r16+r18)*0.5;
				f1[id*np+5] = (r4-r5-2.*r6)*0.5*cs2-(-r3+r12+r14+r17+r18)*0.5;
				f1[id*np+6] = (r4-r5-2.*r6)*0.5*cs2-(r3 -r12-r14+r17+r18)*0.5;
				f1[id*np+7] = (r7+r10+r11+r16)*0.25;
				f1[id*np+8] = (r7+r16-r10-r11)*0.25;
				f1[id*np+9] = (r11+r16-r7-r10)*0.25;
				f1[id*np+10] = (r10+r16-r7-r11)*0.25;
				f1[id*np+11] = (r8+r12+r13+r17)*0.25;
				f1[id*np+12] = (r8+r17-r12-r13)*0.25;
				f1[id*np+13] = (r13+r17-r8-r12)*0.25;
				f1[id*np+14] = (r12+r17-r8-r13)*0.25;
				f1[id*np+15] = (r9+r14+r15+r18)*0.25;
				f1[id*np+16] = (r9+r18-r14-r15)*0.25;
				f1[id*np+17] = (r15+r18-r9-r14)*0.25;
				f1[id*np+18] = (r14+r18-r9-r15)*0.25;

				for(int k=0; k<np; k++)
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
					f2[idn*np+k] = f1[id*np+k];
				}
				kinetic_energy += R*(U*U+V*V+W*W);
				if(isnan(U) || isinf(U))
					check = 1;
			}
	return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  system("mkdir vtk_fluid");
  FILE *data = fopen("data.txt","wt");
	initial_state();
	int check_mach = 0;
	clock_t c_start = clock();
	for(int i=0; i<nsteps; i++)
  {
    check_mach = algo_LB();
    f1 = f2;
    if(check_mach==1)
      goto labelA;
		 if(plot_vtk==true && i%n_out==0)
		 	write_fluid_vtk(i);
		fprintf(data,"%lf    %e\n", (double)i/T_ref, kinetic_energy/kinetic_energy_0);
   	 if(i%1==0)
       printf("Time %lf of %lf. Energy=%e\n", (double)i/T_ref, (double)nsteps/T_ref, kinetic_energy/kinetic_energy_0 );
  }
  labelA:
	clock_t c_end = clock();
	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
