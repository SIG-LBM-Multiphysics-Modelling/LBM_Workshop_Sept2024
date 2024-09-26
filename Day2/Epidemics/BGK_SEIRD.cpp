///----------------------------------------------------------------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
///----------------------------------------------------------------------------------------------------------------------------------
using namespace std;
const bool plot_vtk = false;
const double phys_Lx = 200*1000, //length of the domain in meters
						 phys_Ly = 200*1000, //length of the domain in meters
						 phys_oneday = 1.*24*60*60, // one day in seconds
						 phys_Tmax = 100*phys_oneday, // days to be simulated
						 phys_gamma = pow(5*phys_oneday,-1), // gamma in days^(-1)
						 phys_alpha = pow(7*phys_oneday,-1), // incubation rate in days^(-1)
						 phys_pop = 5E6,
						 phys_rho = phys_pop/(phys_Lx/1000*phys_Ly/1000),
						 phys_dS = 0.0435*(1000*1000)/(24*60*60), // %0.0435 km2/day
						 phys_dE = 0.0198*(1000*1000)/(24*60*60), //
						 phys_dI = 1E-4*(1000*1000)/(24*60*60), //
						 phys_dR = 0.0198*(1000*1000)/(24*60*60), //
						 phys_phi = 0.3;
const int nx = 100, ny = nx, np = 5;
const double cs2 = 1./3., Scale_length = phys_Lx/(nx), oneday = 100, Scale_time = 864, Tmax = phys_Tmax/Scale_time, Scale_gamma = 1./Scale_time, rho0 = phys_pop/((nx)*(ny)),
						 Scale_diffusivity = pow(Scale_length,2)/Scale_time;
const int nsteps = (int)(Tmax)+1, n_out = (int)(Tmax/100);
const vector<int> cx = {0, 1, 0, -1,  0},
									cy = {0, 0, 1,  0, -1},
								 opp = {0, 2, 1,  4,  3};
const vector<double> wf = {1/3., 1/6., 1/6., 1/6., 1/6.};
const double R0 = 3., gamma_ = phys_gamma/Scale_gamma, beta_ = R0*gamma_, alpha = phys_alpha/Scale_gamma, phi = phys_phi;
const double dS = phys_dS/Scale_diffusivity, tauS = dS/cs2+0.5, omegaS = 1./tauS, omega1S = 1.-omegaS;
const double dE = phys_dE/Scale_diffusivity, tauE = dE/cs2+0.5, omegaE = 1./tauE, omega1E = 1.-omegaE;
const double dI = phys_dI/Scale_diffusivity, tauI = dI/cs2+0.5, omegaI = 1./tauI, omega1I = 1.-omegaI;
const double dR = phys_dR/Scale_diffusivity, tauR = dR/cs2+0.5, omegaR = 1./tauR, omega1R = 1.-omegaR;
vector<double> fS1(nx*ny*np, 0.), fS2(nx*ny*np, 0.), rhoS(nx*ny, 0.);
vector<double> fE1(nx*ny*np, 0.), fE2(nx*ny*np, 0.), rhoE(nx*ny, 0.);
vector<double> fI1(nx*ny*np, 0.), fI2(nx*ny*np, 0.), rhoI(nx*ny, 0.);
vector<double> fR1(nx*ny*np, 0.), fR2(nx*ny*np, 0.), rhoR(nx*ny, 0.);
vector<double> rhoD(nx*ny, 0.);
int newx, newy, id, idn, check, idnp;
double RR, RS, RI, RE, RD, feq, N, ftemp, susc, exposed, infected, recovered, dead;
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

	output_file << "SCALARS susceptbiles double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(rhoS[X*ny+Y]>1E-10)
				output_file << rhoS[X*ny+Y]<< "\n";
			else
				output_file << 0 << "\n";
		}

	output_file << "SCALARS exposed double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(rhoE[X*ny+Y]>1E-10)
				output_file << rhoE[X*ny+Y]<< "\n";
			else
				output_file << 0 << "\n";
		}

	output_file << "SCALARS infected double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(rhoI[X*ny+Y]>1E-10)
				output_file << rhoI[X*ny+Y]<< "\n";
			else
				output_file << 0 << "\n";
		}

	output_file << "SCALARS recovered double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(rhoR[X*ny+Y]>1E-10)
				output_file << rhoR[X*ny+Y]<< "\n";
			else
				output_file << 0 << "\n";
		}

	/// Close file
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	srand(time(NULL));
	susc = exposed = infected = recovered = dead = N = 0.;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
    {
			id = x*ny+y;
			rhoI[id] = 0;
			rhoE[id] = rho0*0.01;//*((double)rand()/(RAND_MAX)+0.5);
      rhoS[id] = rho0-rhoE[id];
			rhoR[id] = 0.;
			rhoD[id] = 0.;
			for(int k=0; k<np; k++)
			{
        fS1[id*np+k] = fS2[id*np+k] = wf[k]*rhoS[id];
				fE1[id*np+k] = fE2[id*np+k] = wf[k]*rhoE[id];
				fI1[id*np+k] = fI2[id*np+k] = wf[k]*rhoI[id];
				fR1[id*np+k] = fR2[id*np+k] = wf[k]*rhoR[id];
			}
			susc += rhoS[id];
			exposed += rhoE[id];
			infected += rhoI[id];
			recovered += rhoR[id];
			dead += rhoD[id];
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
int algorithm_lattice_boltzmann()
{
	check = 0;
	susc = exposed = infected = recovered = dead = N = 0.;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			idnp = id*np;
			RS = fS1[idnp+0] + fS1[idnp+1] + fS1[idnp+2] + fS1[idnp+3] + fS1[idnp+4];
			RE = fE1[idnp+0] + fE1[idnp+1] + fE1[idnp+2] + fE1[idnp+3] + fE1[idnp+4];
			RI = fI1[idnp+0] + fI1[idnp+1] + fI1[idnp+2] + fI1[idnp+3] + fI1[idnp+4];
			RR = fR1[idnp+0] + fR1[idnp+1] + fR1[idnp+2] + fR1[idnp+3] + fR1[idnp+4];
			rhoS[id] = RS;
			rhoE[id] = RE;
			rhoI[id] = RI;
			rhoR[id] = RR;
			RD = rhoD[id];
			susc += RS;
			exposed += RE;
			infected += RI;
			recovered += RR;
			dead += RD;
			N = RS+RE+RI+RR+RD;
			for(int k=0; k<np; k++)
      {
      	newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
				idn = newx*ny+newy;
				fS2[idn*np+k] = omega1S*fS1[idnp+k] + omegaS*wf[k]*RS - wf[k]*beta_*RS*RI/N;
				fE2[idn*np+k] = omega1E*fE1[idnp+k] + omegaE*wf[k]*RE + wf[k]*(beta_*RS*RI/N-alpha*RE);
				fI2[idn*np+k] = omega1I*fI1[idnp+k] + omegaI*wf[k]*RI + wf[k]*(alpha*RE-gamma_*RI);
				fR2[idn*np+k] = omega1R*fR1[idnp+k] + omegaR*wf[k]*RR + wf[k]*gamma_*(1.-phi)*RI;
			}
			rhoD[id] += gamma_*phi*RI;
			if(isnan(RS))
				check = 0;
		}
    return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	FILE *data = fopen("dataBGK.txt","wt");
	system("mkdir vtk_fluid");
	int check_mach = 0, t;
	initial_state();
	printf("%e %e\n", Scale_length, Scale_time);
	clock_t c_start = clock();
	for(t=0; t<nsteps; t++)
  {
    check_mach = algorithm_lattice_boltzmann();
   	//boundary();
		fS1 = fS2;
		fE1 = fE2;
		fI1 = fI2;
		fR1 = fR2;
		if(plot_vtk==true && t%n_out==0)
			write_fluid_vtk(t);
		if(t%n_out==0)
			printf("Day %lf of %lf. Infected=%e\n", t*Scale_time/phys_oneday, nsteps*Scale_time/phys_oneday, infected/(susc+exposed+infected+recovered+dead));
		fprintf(data,"%lf    %e    %e    %e    %e    %e    %e\n", (double)t*Scale_time/phys_oneday, (susc+exposed+infected+recovered+dead)/(susc+exposed+infected+recovered+dead), susc/(susc+exposed+infected+recovered+dead), exposed/(susc+exposed+infected+recovered+dead), infected/(susc+exposed+infected+recovered+dead), recovered/(susc+exposed+infected+recovered+dead), dead/(susc+exposed+infected+recovered+dead));
		if(check_mach==1)
      goto labelA;
  }
  labelA:
  clock_t c_end = clock();
	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	printf("Time = %lf\n", time_elapsed_ms);
	fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
