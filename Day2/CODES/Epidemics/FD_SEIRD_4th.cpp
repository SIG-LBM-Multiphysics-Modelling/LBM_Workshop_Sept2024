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
const int nx = 100, ny = nx, np = 9;
const double cs2 = 1./3., Scale_length = phys_Lx/(nx), oneday = 100, Scale_time = 864, Tmax = phys_Tmax/Scale_time, Scale_gamma = 1./Scale_time, rho0 = phys_pop/((nx)*(ny)),
						 Scale_diffusivity = pow(Scale_length,2)/Scale_time;
const int nsteps = (int)(Tmax)+1, n_out = (int)(Tmax/100);
const vector<int> cx = {0, 1, 0, -1,  0, 1, -1, -1, 1},
									cy = {0, 0, 1,  0, -1, 1, 1, -1, -1};
const vector<double> wf = {4/9., 1/9., 1/9., 1/9., 1/9., 1/36., 1/36., 1/36., 1/36.};
const double R0 = 3., phys_beta = R0*phys_gamma;
vector<double> rhoS(nx*ny, 0.), rhoS_old(nx*ny, 0.), rhoE(nx*ny, 0.), rhoE_old(nx*ny, 0.), rhoI(nx*ny, 0.), rhoI_old(nx*ny, 0.), rhoR(nx*ny, 0.), rhoR_old(nx*ny, 0.), rhoD(nx*ny, 0.), rhoD_old(nx*ny, 0.);
int newx, newy, id, idn, check, idnp;
double laplR, laplS, laplI, laplE, RD, feq, N, ftemp, susc, exposed, infected, recovered, dead;
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
			rhoE[id] = phys_rho*0.01;//*((double)rand()/(RAND_MAX)+0.5);
      rhoS[id] = phys_rho-rhoE[id];
			rhoR[id] = 0.;
			rhoD[id] = 0.;

			susc += rhoS[id];
			exposed += rhoE[id];
			infected += rhoI[id];
			recovered += rhoR[id];
			dead += rhoD[id];
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
int FD()
{
	check = 0;
	rhoS_old = rhoS;
	rhoE_old = rhoE;
	rhoI_old = rhoI;
	rhoR_old = rhoR;
	rhoD_old = rhoD;
	susc = exposed = infected = recovered = dead = N = 0.;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;;
			laplS = laplE = laplI = laplR = 0.;
			for(int k=1; k<np; k++)
      {
      	newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
				idn = newx*ny+newy;
        laplS += 6.*wf[k]*(rhoS_old[idn]-rhoS_old[id]);
				laplE += 6.*wf[k]*(rhoE_old[idn]-rhoE_old[id]);
				laplI += 6.*wf[k]*(rhoI_old[idn]-rhoI_old[id]);
				laplR += 6.*wf[k]*(rhoR_old[idn]-rhoR_old[id]);
      }
      N = rhoS_old[id]+rhoE_old[id]+rhoI_old[id]+rhoR_old[id]+rhoD_old[id];
			rhoS[id] = rhoS_old[id] + 0.25*Scale_time*(-phys_beta*rhoS_old[id]*rhoI_old[id]/N + phys_dS*laplS/pow(Scale_length,2));
			rhoE[id] = rhoE_old[id] + 0.25*Scale_time*(phys_beta*rhoS_old[id]*rhoI_old[id]/N - phys_alpha*rhoE_old[id] + phys_dE*laplE/pow(Scale_length,2));
			rhoI[id] = rhoI_old[id] + 0.25*Scale_time*(phys_alpha*rhoE_old[id] - phys_gamma*rhoI_old[id] + phys_dI*laplI/pow(Scale_length,2));
			rhoR[id] = rhoR_old[id] + 0.25*Scale_time*(phys_gamma*(1.-phys_phi)*rhoI_old[id] + phys_dR*laplR/pow(Scale_length,2));
			rhoD[id] = rhoD_old[id] + 0.25*Scale_time*phys_gamma*phys_phi*rhoI_old[id];
		}
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;;
			laplS = laplE = laplI = laplR = 0.;
			for(int k=1; k<np; k++)
	    {
	    	newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
				idn = newx*ny+newy;
	      laplS += 6.*wf[k]*(rhoS[idn]-rhoS[id]);
				laplE += 6.*wf[k]*(rhoE[idn]-rhoE[id]);
				laplI += 6.*wf[k]*(rhoI[idn]-rhoI[id]);
				laplR += 6.*wf[k]*(rhoR[idn]-rhoR[id]);
	    }
	    N = rhoS[id]+rhoE[id]+rhoI[id]+rhoR[id]+rhoD[id];
			rhoS_old[id] = rhoS[id] + 0.25*Scale_time*(-phys_beta*rhoS[id]*rhoI[id]/N + phys_dS*laplS/pow(Scale_length,2));
			rhoE_old[id] = rhoE[id] + 0.25*Scale_time*(phys_beta*rhoS[id]*rhoI[id]/N - phys_alpha*rhoE[id] + phys_dE*laplE/pow(Scale_length,2));
			rhoI_old[id] = rhoI[id] + 0.25*Scale_time*(phys_alpha*rhoE[id] - phys_gamma*rhoI[id] + phys_dI*laplI/pow(Scale_length,2));
			rhoR_old[id] = rhoR[id] + 0.25*Scale_time*(phys_gamma*(1.-phys_phi)*rhoI[id] + phys_dR*laplR/pow(Scale_length,2));
			rhoD_old[id] = rhoD[id] + 0.25*Scale_time*phys_gamma*phys_phi*rhoI[id];
		}

	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;;
			laplS = laplE = laplI = laplR = 0.;
			for(int k=1; k<np; k++)
      {
      	newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
				idn = newx*ny+newy;
        laplS += 6.*wf[k]*(rhoS_old[idn]-rhoS_old[id]);
				laplE += 6.*wf[k]*(rhoE_old[idn]-rhoE_old[id]);
				laplI += 6.*wf[k]*(rhoI_old[idn]-rhoI_old[id]);
				laplR += 6.*wf[k]*(rhoR_old[idn]-rhoR_old[id]);
      }
      N = rhoS_old[id]+rhoE_old[id]+rhoI_old[id]+rhoR_old[id]+rhoD_old[id];
			rhoS[id] = rhoS_old[id] + 0.25*Scale_time*(-phys_beta*rhoS_old[id]*rhoI_old[id]/N + phys_dS*laplS/pow(Scale_length,2));
			rhoE[id] = rhoE_old[id] + 0.25*Scale_time*(phys_beta*rhoS_old[id]*rhoI_old[id]/N - phys_alpha*rhoE_old[id] + phys_dE*laplE/pow(Scale_length,2));
			rhoI[id] = rhoI_old[id] + 0.25*Scale_time*(phys_alpha*rhoE_old[id] - phys_gamma*rhoI_old[id] + phys_dI*laplI/pow(Scale_length,2));
			rhoR[id] = rhoR_old[id] + 0.25*Scale_time*(phys_gamma*(1.-phys_phi)*rhoI_old[id] + phys_dR*laplR/pow(Scale_length,2));
			rhoD[id] = rhoD_old[id] + 0.25*Scale_time*phys_gamma*phys_phi*rhoI_old[id];
		}
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;;
			laplS = laplE = laplI = laplR = 0.;
			for(int k=1; k<np; k++)
	    {
	    	newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
				idn = newx*ny+newy;
	      laplS += 6.*wf[k]*(rhoS[idn]-rhoS[id]);
				laplE += 6.*wf[k]*(rhoE[idn]-rhoE[id]);
				laplI += 6.*wf[k]*(rhoI[idn]-rhoI[id]);
				laplR += 6.*wf[k]*(rhoR[idn]-rhoR[id]);
	    }
	    N = rhoS[id]+rhoE[id]+rhoI[id]+rhoR[id]+rhoD[id];
			rhoS_old[id] = rhoS[id] + 0.25*Scale_time*(-phys_beta*rhoS[id]*rhoI[id]/N + phys_dS*laplS/pow(Scale_length,2));
			rhoE_old[id] = rhoE[id] + 0.25*Scale_time*(phys_beta*rhoS[id]*rhoI[id]/N - phys_alpha*rhoE[id] + phys_dE*laplE/pow(Scale_length,2));
			rhoI_old[id] = rhoI[id] + 0.25*Scale_time*(phys_alpha*rhoE[id] - phys_gamma*rhoI[id] + phys_dI*laplI/pow(Scale_length,2));
			rhoR_old[id] = rhoR[id] + 0.25*Scale_time*(phys_gamma*(1.-phys_phi)*rhoI[id] + phys_dR*laplR/pow(Scale_length,2));
			rhoD_old[id] = rhoD[id] + 0.25*Scale_time*phys_gamma*phys_phi*rhoI[id];

			susc += rhoS_old[id];
			exposed += rhoE_old[id];
			infected += rhoI_old[id];
			recovered += rhoR_old[id];
			dead += rhoD_old[id];

			if(isnan(rhoS_old[id]))
				check = 0;
		}

	rhoS = rhoS_old; 
	rhoE = rhoE_old; 
  rhoI = rhoI_old; 
  rhoR = rhoR_old; 
  rhoD = rhoD_old; 

  return check;
}

///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	FILE *data = fopen("dataFD.txt","wt");
	system("mkdir vtk_fluid");
	int check_mach = 0, t;
	initial_state();
	printf("%e %e\n", Scale_length, Scale_time);
	clock_t c_start = clock();
	for(t=0; t<nsteps; t++)
  {
    check_mach = FD();
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
	float time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";
	fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
