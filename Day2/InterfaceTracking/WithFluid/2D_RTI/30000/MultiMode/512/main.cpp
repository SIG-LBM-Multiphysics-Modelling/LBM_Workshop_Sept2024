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
const int nx = 64, ny = 2*nx, nz = 1, np = 19;
//          						0  1   2  3   4  5   6  7   8   9  10 11  12  13  14 15  16  17  18
vector<const int> cx = {0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1, 0,  0,  0,  0},
									cy = {0, 0,  0, 1, -1, 0,  0, 1, -1, -1,  1, 0,  0,  0,  0, 1, -1,  1, -1},
									cz = {0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1, -1, -1,  1, 1, -1, -1,  1},
								 opp = {0, 2,  1, 4,  3, 6,  5, 8,  7, 10,  9, 12, 11, 14, 13, 16, 15, 18, 17};
const double Reynolds = 3000., At = 0.1, Ca = 0.26, Pe = 500.;
const double cs2 = 1./3., cs4 = cs2*cs2, xi = 5., d = (double)nx, U_ref = 0.04, gravity = pow(U_ref,2)/d, nu = sqrt(d*gravity)*d/Reynolds, T_ref = sqrt(d/gravity/At);
vector<const double> wf = {1/3., 1/18., 1/18., 1/18., 1/18., 1/18., 1/18., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36.};
const int nsteps = (int)(3.*T_ref+1), n_out = (int)(0.5*T_ref);
const double rhoL = 1., rhoH = rhoL*(1.+At)/(1.-At), niH = nu, tauH = niH/cs2, niL = nu, miL = niL*rhoL, miH = niH*rhoH, tauL = niL/cs2, sigma = niH*U_ref/Ca;
vector<double> f1(nx*ny*nz*np,0.), f2(nx*ny*nz*np,0.), rho(nx*ny*nz,0.), press(nx*ny*nz,0.), u(nx*ny*nz,0.), v(nx*ny*nz,0.), w(nx*ny*nz,0.), temp_pop(np,0.);
double feq, tau, omega, omega1, ni, Fpx, Fpy, Fpz, Fsx, Fsy, Fsz, Fmx, Fmy, Fmz, mu, Press, Press0, U0, V0, W0, mi;
double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18;
double r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18;
double R, U, V, W, ftemp, CX, CY, CZ, A, B, C, U2, V2, W2, Fx, Fy, Fz, U3, V3, W3, U4, V4, W4, P1;
int newx, newy, newz, id, idn;
// PHASE
double k1_g, k2_g, k3_g, k4_g, k5_g, k6_g, k7_g, k8_g, k9_g, k10_g, k11_g, k12_g, k13_g, k14_g, k15_g, k16_g, k17_g, k18_g;
const double PhiH = 1., PhiL = 0., Phi0 = 0.5*(PhiH+PhiL), beta = 12.*sigma/xi, kappa = 3.*sigma*xi/2.;
const double M = U_ref*d/Pe, tau_phase = M/cs2+0.5, omega_phase = 1./tau_phase, omega_phase1 = 1.-omega_phase;
double grad_phix, grad_phiy, grad_phiz, laplPhi, Phi, gtemp, P, Fx_phase, Fy_phase, Fz_phase, Nx, Ny, Nz, Phii, thisPhi, Phi_prev;
vector<double> g1(nx*ny*nz*np,0.), g2(nx*ny*nz*np,0.), phase(nx*ny*nz,0.), phase_old(nx*ny*nz,0.), temp_pop_phase(np,0.);
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
	// output_file << "SCALARS density float 1\n";
	// output_file << "LOOKUP_TABLE default\n";
	// for(int Z = 0; Z < nz; ++Z)
	// 	for(int Y = 0; Y < ny ; ++Y)
	// 		for(int X = 0; X < nx; ++X)
	// 		{
	// 				id = (X*ny+Y)*nz+Z;
	// 				if(rho[id]>1E-12)
	// 					output_file << rho[id] << "\n";
	// 				else
	// 					output_file << "0 \n";
	// 		}
	//
	// output_file << "VECTORS velocity_vector float\n";
	// for(int Z = 0; Z < nz; ++Z)
	// 	for(int Y = 0; Y < ny ; ++Y)
	// 		for(int X = 0; X < nx; ++X)
	// 		{
	// 			id = (X*ny+Y)*nz+Z;
	// 			output_file << u[id] << " " << v[id] << " " << w[id] << "\n";
	// 		}

	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	double extra_term, X, Z;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
    	for(int z=0; z<nz; z++)
			{
				X = double(x) / (2*double(nx-1));
				Z = double(z) / double(nz-1);
				id = (x*ny+y)*nz+z;
				if(y>2.*nx+0.1*nx*(cos(2.*M_PI*X)))
					phase[id] = PhiH;
			  else
					phase[id] = PhiL;
			  rho[id] = rhoL+(phase[id]-PhiL)/(PhiH-PhiL)*(rhoH-rhoL);
				U = u[id] = 0.;
      	V = v[id] = 0.;
	      W = w[id] = 0.;
				C = -1.5*(U*U+V*V+W*W);
				Press = press[id] = 0.;
				for(int k=0; k<np; k++)
				{
	        A = U*cx[k]+V*cy[k]+W*cz[k];
	        B = 4.5*A*A;
	        f1[id*np+k]= wf[k]*(Press+3.*A+B+C);
					g1[id*np+k]= wf[k]*phase[id]*(1.+3.*A+B+C);
				}
			}
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

	if(y==0)
	{
		grad_phiy = phase_old[(x*ny+y+1)*nz+z]-phase_old[id];
		grad_phix = -0.5*phase_old[(((x-1+nx)%nx)*ny+y)*nz+z] + 0.5*phase_old[(((x+1+nx)%nx)*ny+y)*nz+z];
		grad_phiz = -0.5*phase_old[(x*ny+y)*nz+(z-1+nz)%nz] + 0.5*phase_old[(x*ny+y)*nz+(z+1+nz)%nz];
	}
	if(y==ny-1)
	{
		grad_phiy = -phase_old[(x*ny+y-1)*nz+z]+phase_old[id];
		grad_phix = -0.5*phase_old[(((x-1+nx)%nx)*ny+y)*nz+z] + 0.5*phase_old[(((x+1+nx)%nx)*ny+y)*nz+z];
		grad_phiz = -0.5*phase_old[(x*ny+y)*nz+(z-1+nz)%nz] + 0.5*phase_old[(x*ny+y)*nz+(z+1+nz)%nz];
	}
}
///----------------------------------------------------------------------------------------------------------------------------------
int algo_LB()
{
	int hh = 0;
	phase_old = phase;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
    	for(int z=0; z<nz; z++)
			{
				id = (x*ny+y)*nz+z;
				compute_grad_phi(x, y, z);
				U0 = u[id];
				V0 = v[id];
				W0 = w[id];
				Press0 = press[id];
				Phi_prev = phase_old[id];
				U = V = W = Press = Phi = 0.;
				Fmx = Fmy = Fmz = 0.;
				for(int k=0; k<np; ++k)
				{
					temp_pop[k] = ftemp = f1[id*np+k];
	        Press += ftemp;
	        U += ftemp*cx[k];
	        V += ftemp*cy[k];
					W += ftemp*cz[k];

					temp_pop_phase[k] = gtemp = g1[id*np+k];
					gtemp = g1[id*np+k];
					Phi += gtemp;

	        A = U0*cx[k]+V0*cy[k]+W0*cz[k];
	        feq = wf[k]*(Press0+3.*A+4.5*A*A-1.5*(U0*U0+V0*V0+W0*W0));
	        Fmx += cx[k]*cy[k]*(ftemp-feq)*grad_phiy+cx[k]*cz[k]*(ftemp-feq)*grad_phiz;
	        Fmy += cy[k]*cx[k]*(ftemp-feq)*grad_phix+cy[k]*cz[k]*(ftemp-feq)*grad_phiz;
					Fmz += cz[k]*cx[k]*(ftemp-feq)*grad_phix+cz[k]*cy[k]*(ftemp-feq)*grad_phiy;
				}
				phase[id] = Phi;
				press[id] = Press;
				Fpx = -Press0*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix;
				Fpy = -Press0*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy;
				Fpz = -Press0*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiz;
				mu = 4.*beta*(Phi_prev-PhiL)*(Phi_prev-PhiH)*(Phi_prev-Phi0)-kappa*laplPhi;
				Fsx = mu*grad_phix;
				Fsy = mu*grad_phiy;
				Fsz = mu*grad_phiz;
				mi = miL+(Phi_prev-PhiL)/(PhiH-PhiL)*(miH-miL);
				rho[id] = R = rhoL+(Phi_prev-PhiL)/(PhiH-PhiL)*(rhoH-rhoL);
				tau = mi/R/cs2;
				omega = 1./(tau+0.5);
				omega1 = 1.-omega;
				Fmx *= -ni/((tau+0.5)*cs2)*(rhoH-rhoL)/(PhiH-PhiL);
				Fmy *= -ni/((tau+0.5)*cs2)*(rhoH-rhoL)/(PhiH-PhiL);
				Fmz *= -ni/((tau+0.5)*cs2)*(rhoH-rhoL)/(PhiH-PhiL);

				Fx = Fpx+Fsx+Fmx;
				Fy = Fpy+Fsy+Fmy;
				Fz = Fpz+Fsz+Fmz;
				Fy += -(R-0.5*(rhoH+rhoL))*gravity;
				//Fy += -R*gravity;
				U += 0.5*Fx/R;
				V += 0.5*Fy/R;
				W += 0.5*Fz/R;
				u[id] = U;
				v[id] = V;
				w[id] = W;
				if(y==0 || y==ny-1)
					u[id] = v[id] = U = V = w[id] = W = 0.;
				Nx  = grad_phix/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2)+pow(grad_phiz,2))+1E-12);
				Ny  = grad_phiy/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2)+pow(grad_phiz,2))+1E-12);
				Nz  = grad_phiz/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2)+pow(grad_phiz,2))+1E-12);
				U2 = U*U;
				V2 = V*V;
				W2 = W*W;
				if(fabs(U)>1.)
          hh = 1;
        ///compute moments "k"
				/*k5 = k6 = k7 = k8 = k9 = 0.;
				k1_g = k2_g = k3_g = 0.;
        for(int k=0; k<np; ++k)
        {
          CX = cx[k]-U;
          CY = cy[k]-V;
          CZ = cz[k]-W;
					ftemp = temp_pop[k];
					//k5 += ftemp*(CX*CX-CY*CY);
					//k6 += ftemp*(CY*CY-CZ*CZ);
					//k7 += ftemp*CX*CY;
					//k8 += ftemp*CX*CZ;
					//k9 += ftemp*CY*CZ;

					gtemp = temp_pop_phase[k];
					//k1_g += ftemp*CX;
					//k2_g += ftemp*CY;
					//k3_g += ftemp*CZ;
				}*/
				///collide moments
				/*r5 = temp_pop[1]+temp_pop[2]-temp_pop[3]-temp_pop[4]+temp_pop[11]+temp_pop[12]+temp_pop[13]+temp_pop[14]-temp_pop[15]-temp_pop[16]-temp_pop[17]-temp_pop[18];
				r6 = temp_pop[3]+temp_pop[4]-temp_pop[5]-temp_pop[6]+temp_pop[7]+temp_pop[8]+temp_pop[9]+temp_pop[10]-temp_pop[11]-temp_pop[12]-temp_pop[13]-temp_pop[14];
				r7 = temp_pop[7]+temp_pop[8]-temp_pop[9]-temp_pop[10];
				r8 = temp_pop[11]+temp_pop[12]-temp_pop[13]-temp_pop[14];
				r9 = temp_pop[15]+temp_pop[16]-temp_pop[17]-temp_pop[18];
				k5 = r5+(U2-V2)*(Press-2.);
				k6 = r6+(V2-W2)*(Press-2.);
				k7 = r7+U*V*(Press-2.);
				k8 = r8+U*W*(Press-2.);
				k9 = r9+V*W*(Press-2.);
				P1 = Press-1.;
				k1 = Fx*0.5-U*P1;
				k2 = Fy*0.5-V*P1;
				k3 = Fz*0.5-W*P1;
				k4 = Press*(1.+U2+V2+W2)-U2-V2-W2;
				k5 = omega1*k5+omega*(U2-V2)*P1;
				k6 = omega1*k6+omega*(V2-W2)*P1;
				k7 = omega1*k7+U*V*omega*P1;
				k8 = omega1*k8+U*W*omega*P1;
				k9 = omega1*k9+V*W*omega*P1;
				k10 = Fy*0.5*cs2-V*(3.*U2+1.)*P1*cs2;
				k11 = Fx*0.5*cs2-U*(3.*V2+1.)*P1*cs2;
				k12 = Fz*0.5*cs2-W*(3.*U2+1.)*P1*cs2;
				k13 = Fx*0.5*cs2-U*(3.*W2+1.)*P1*cs2;
				k14 = Fz*0.5*cs2-W*(3.*V2+1.)*P1*cs2;
				k15 = Fy*0.5*cs2-V*(3.*W2+1.)*P1*cs2;
				k16 = Press*(cs4+cs2*(U2+V2)+U2*V2)-cs2*(U2+V2)-U2*V2;
				k17 = Press*(cs4+cs2*(U2+W2)+U2*W2)-cs2*(U2+W2)-U2*W2;
				k18 = Press*(cs4+cs2*(V2+W2)+V2*W2)-cs2*(V2+W2)-V2*W2;
				r0 = Press;
				r1 = k1+Press*U;
				r2 = k2+Press*V;
				r3 = k3+Press*W;
				r4 = Press*(U2+V2+W2)+2.*k1*U+2.*k2*V+2.*k3*W+k4;
				r5 = Press*(U2-V2)+2.*k1*U-2.*k2*V+k5;
				r6 = Press*(V2-W2)+2.*k2*V-2.*k3*W+k6;
				r7 = k7+U*k2+V*k1+Press*U*V;
				r8 = k8+U*k3+W*k1+Press*U*W;
				r9 = k9+V*k3+W*k2+Press*V*W;
				r10 = k10+2.*U*k7+V*k4*cs2+2.*V*k5*cs2+V*k6*cs2+U2*k2+2.*U*V*k1+Press*U2*V;
				r11 = k11+U*k4*cs2-U*k5*cs2+U*k6*cs2+2.*V*k7+V2*k1+2.*U*V*k2+Press*U*V2;
				r12 = k12+2.*U*k8+W*k4*cs2+2.*W*k5*cs2+W*k6*cs2+U2*k3+2.*U*W*k1+Press*U2*W;
				r13 = k13+U*k4*cs2-U*k5*cs2-2.*U*k6*cs2+2.*W*k8+W2*k1+2.*U*W*k3+Press*U*W2;
				r14 = k14+2.*V*k9+W*k4*cs2-W*k5*cs2+W*k6*cs2+V2*k3+2.*V*W*k2+Press*V2*W;
				r15 = k15+V*k4*cs2-V*k5*cs2-2.*V*k6*cs2+2.*W*k9+W2*k2+2.*V*W*k3+Press*V*W2;
				r16 = k16+2.*U*k11+2.*V*k10+U2*k4*cs2-U2*k5*cs2+U2*k6*cs2+V2*k4*cs2+2.*V2*k5*cs2+V2*k6*cs2+Press*U2*V2+4.*U*V*k7+2.*U*V2*k1+2.*U2*V*k2;
				r17 = k17+2.*U*k13+2.*W*k12+U2*k4*cs2-U2*k5*cs2-2.*U2*k6*cs2+W2*k4*cs2+2.*W2*k5*cs2+W2*k6*cs2+Press*U2*W2+4.*U*W*k8+2.*U*W2*k1+2.*U2*W*k3;
				r18 = k18+2.*V*k15+2.*W*k14+V2*k4*cs2-V2*k5*cs2-2.*V2*k6*cs2+W2*k4*cs2-W2*k5*cs2+W2*k6*cs2+Press*V2*W2+4.*V*W*k9+2.*V*W2*k2+2.*V2*W*k3;
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
				f1[id*np+18] = (r14+r18-r9-r15)*0.25;*/

				Fx_phase = cs2*(1.-4.*pow((Phi-Phi0),2))/xi*Nx;
				Fy_phase = cs2*(1.-4.*pow((Phi-Phi0),2))/xi*Ny;
				Fz_phase = cs2*(1.-4.*pow((Phi-Phi0),2))/xi*Nz;
				r1 = temp_pop_phase[1]-temp_pop_phase[2]+temp_pop_phase[7]-temp_pop_phase[8]+temp_pop_phase[9]-temp_pop_phase[10]+temp_pop_phase[11]-temp_pop_phase[12]+temp_pop_phase[13]-temp_pop_phase[14];
				r2 = temp_pop_phase[3]-temp_pop_phase[4]+temp_pop_phase[7]-temp_pop_phase[8]-temp_pop_phase[9]+temp_pop_phase[10]+temp_pop_phase[15]-temp_pop_phase[16]+temp_pop_phase[17]-temp_pop_phase[18];
				r3 = temp_pop_phase[5]-temp_pop_phase[6]+temp_pop_phase[11]-temp_pop_phase[12]-temp_pop_phase[13]+temp_pop_phase[14]+temp_pop_phase[15]-temp_pop_phase[16]-temp_pop_phase[17]+temp_pop_phase[18];
				k1_g = r1-Phi*U;
				k2_g = r2-Phi*V;
				k3_g = r3-Phi*W;
				k1_g = k1_g*omega_phase1+(1.-omega_phase*0.5)*Fx_phase;
				k2_g = k2_g*omega_phase1+(1.-omega_phase*0.5)*Fy_phase;
				k3_g = k3_g*omega_phase1+(1.-omega_phase*0.5)*Fz_phase;
				k4_g = Phi;
				k5_g = 0.;
				k6_g = 0.;
				k7_g = 0.;
				k8_g = 0.;
				k9_g = 0.;
				k10_g = Fy_phase*0.5*cs2;
				k11_g = Fx_phase*0.5*cs2;
				k12_g = Fz_phase*0.5*cs2;
				k13_g = Fx_phase*0.5*cs2;
				k14_g = Fz_phase*0.5*cs2;
				k15_g = Fy_phase*0.5*cs2;
				k16_g = Phi*cs4;
				k17_g = Phi*cs4;
				k18_g = Phi*cs4;
				r0 = Phi;
				r1 = k1_g+Phi*U;
				r2 = k2_g+Phi*V;
				r3 = k3_g+Phi*W;
				r4 = Phi*(U2+V2+W2)+2.*k1_g*U+2.*k2_g*V+2.*k3_g*W+k4_g;
				r5 = Phi*(U2-V2)+2.*k1_g*U-2.*k2_g*V;
				r6 = Phi*(V2-W2)+2.*V*k2_g-2.*W*k3_g;
				r7 = U*k2_g+V*k1_g+Phi*U*V;
				r8 = U*k3_g+W*k1_g+Phi*U*W;
				r9 = V*k3_g+W*k2_g+Phi*V*W;
				r10 = k10_g+V*k4_g*cs2+U2*k2_g+2.*U*V*k1_g+Phi*U2*V;
				r11 = k11_g+U*k4_g*cs2+V2*k1_g+2.*U*V*k2_g+Phi*U*V2;
				r12 = k12_g+W*k4_g*cs2+U2*k3_g+2.*U*W*k1_g+Phi*U2*W;
				r13 = k13_g+U*k4_g*cs2+W2*k1_g+2.*U*W*k3_g+Phi*U*W2;
				r14 = k14_g+W*k4_g*cs2+V2*k3_g+2.*V*W*k2_g+Phi*V2*W;
				r15 = k15_g+V*k4_g*cs2+W2*k2_g+2.*V*W*k3_g+Phi*V*W2;
				r16 = k16_g+2.*U*k11_g+2.*V*k10_g+k4_g*(U2*cs2+V2*cs2)+Phi*U2*V2+2.*U*V2*k1_g+2.*U2*V*k2_g;
				r17 = k17_g+2.*U*k13_g+2.*W*k12_g+k4_g*(U2*cs2+W2*cs2)+Phi*U2*W2+2.*U*W2*k1_g+2.*U2*W*k3_g;
				r18 = k18_g+Phi*V2*W2+2.*k3_g*V2*W+k4_g*V2*cs2+2.*k2_g*V*W2+2.*k15_g*V+k4_g*W2*cs2+2.*k14_g*W;
				g1[id*np+0] = r0-r4+r16+r17+r18;
				g1[id*np+1] = (r4+2.*r5+r6)*0.5*cs2-(-r1+r11+r13+r16+r17)*0.5;
				g1[id*np+2] = (r4+2.*r5+r6)*0.5*cs2+(r11+r13-r1-r16-r17)*0.5;
				g1[id*np+3] = (r4-r5+r6)*0.5*cs2-(-r2+r10+r15+r16+r18)*0.5;
				g1[id*np+4] = (r4-r5+r6)*0.5*cs2-(r2-r10-r15+r16+r18)*0.5;
				g1[id*np+5] = (r4-r5-2.*r6)*0.5*cs2-(-r3+r12+r14+r17+r18)*0.5;
				g1[id*np+6] = (r4-r5-2.*r6)*0.5*cs2-(r3 -r12-r14+r17+r18)*0.5;
				g1[id*np+7] = (r7+r10+r11+r16)*0.25;
				g1[id*np+8] = (r7+r16-r10-r11)*0.25;
				g1[id*np+9] = (r11+r16-r7-r10)*0.25;
				g1[id*np+10] = (r10+r16-r7-r11)*0.25;
				g1[id*np+11] = (r8+r12+r13+r17)*0.25;
				g1[id*np+12] = (r8+r17-r12-r13)*0.25;
				g1[id*np+13] = (r13+r17-r8-r12)*0.25;
				g1[id*np+14] = (r12+r17-r8-r13)*0.25;
				g1[id*np+15] = (r9+r14+r15+r18)*0.25;
				g1[id*np+16] = (r9+r18-r14-r15)*0.25;
				g1[id*np+17] = (r15+r18-r9-r14)*0.25;
				g1[id*np+18] = (r14+r18-r9-r15)*0.25;
				C = -1.5*(U2+V2+W2);
				for(int k=0; k<np; ++k)
				{
					A = U*cx[k]+V*cy[k]+W*cz[k];
	        B = 4.5*A*A;
					f1[id*np+k] *= omega1;
					f1[id*np+k] += omega*wf[k]*(Press+3.*A+B+C);
					f1[id*np+k] += (1.-0.5*omega)*wf[k]*( Fx*(3*(cx[k]-U)+9*A*cx[k])+
					 																			Fy*(3*(cy[k]-V)+9*A*cy[k])+
					 																			Fz*(3*(cz[k]-W)+9*A*cz[k]));
					// g1[id*np+k] *= omega_phase1;
					// g1[id*np+k] += omega_phase*wf[k]*phase[id]*(1.+3.*A+B+C);
					// g1[id*np+k] += (1.-0.5*omega_phase)*wf[k]*( Fx_phase*(3*(cx[k]-U)+9*A*cx[k])+
					// 																						Fy_phase*(3*(cy[k]-V)+9*A*cy[k])+
					// 																						Fz_phase*(3*(cz[k]-W)+9*A*cz[k]));
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
					g2[idn*np+k] = g1[id*np+k];
				}
			}
  return hh;
}
///----------------------------------------------------------------------------------------------------------------------------------
void boundary_conditions()
{
  for(int x=0; x<nx; x++)
		for(int z=0; z<nz; z++)
			for(int k=0; k<np; k++)
			{
				id = (x*ny+0)*nz+z;
				if(cy[k]>0)
				{
					f2[id*np+k] = f1[id*np+opp[k]];
					g2[id*np+k] = g1[id*np+opp[k]];
				}
				id = (x*ny+ny-1)*nz+z;
				if(cy[k]<0)
				{
					f2[id*np+k] = f1[id*np+opp[k]];
					g2[id*np+k] = g1[id*np+opp[k]];
				}
			}
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
    check_mach = algo_LB();
		boundary_conditions();
		f1 = f2;
		g1 = g2;
    //fprintf(data,"%lf    %e    %e\n", (double)i/T_ref, kinetic_energy, enstrophy);
    if(check_mach==1)
      goto labelA;
		if(plot_vtk==true && i%n_out==0)
			write_fluid_vtk(i);
   	if(i%1==0)
      printf("Iteration %lf of %lf. \n", (double)i/T_ref, (double)nsteps/T_ref);
  }
	  labelA:
  clock_t c_end = clock();
	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
