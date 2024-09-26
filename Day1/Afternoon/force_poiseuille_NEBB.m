% This code is designed for the LBM Workshop 2024 @ Universtiy of
% Manchester.
% This code is provided under the MIT license.
% Author: Goncalo Silva
%
% =========================================================================
% Example matlab code for computing a Poiseuille flow with BB boundary


clear all
close all
clc

% simulation parameters
NX=3;             % channel length
NY=16;             % channel width
NSTEPS=1e6;     % number of simulation time steps



rho0=1;                 % rest density state
tau=0.9;                % relaxation time (BGK model)
nu=1/3*(tau-1/2);         % kinematic shear viscosity
Re=10.0;                 % Reynolds number; scaling parameter in simulation (Re=NY*u_max/nu)
u_max=Re*nu/NY;              % maximum velocity

forcex = 8.*u_max*nu./(NY.^2); % body force component x
forcey=0;                        % body force component y

% Lattice parameters; note zero direction is last
NPOP=9;                                         % number of velocities
w  = [1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9]; % weights
cx = [1 0 -1  0 1 -1 -1  1 0];                  % velocities, x components
cy = [0 1  0 -1 1  1 -1 -1 0];                  % velocities, y components


% Node locations
x = (1:NX);
y = (1:NY);

% Analytical solution: Poiseuille velocity
u_analy=-4*u_max/(NY^2).*(y-1).*(y-NY);


% initialize populations
feq=zeros(NX,NY,NPOP);
for k=1:NPOP
    feq(:,:,k)=w(k); % assuming density equal one and zero velocity initial state
end

f=feq;
fcoll=feq;
Source=w.*3.*(cx.*forcex+cy.*forcey);

% Initialize macroscopic populations
rho=ones(NX,NY);
u=zeros(NX,NY);
v=zeros(NX,NY);

% convergence parameters
tol=1e-12;      % tolerance to steady state convergence
teval=100;      % time step to evaluate convergence
u_old=u;

% initalize clock
tstart = tic;

% Main algorithm
for t=1:NSTEPS
    % Compute macroscopic quantities
    
    
    for j=1:NY
        for i=1:NX
            % density
            rho(i,j)=f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+...
                f(i,j,7)+f(i,j,8)+f(i,j,9);
            % velocity components
            u(i,j)=1./rho(i,j).*(f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8));
            v(i,j)=1./rho(i,j).*(f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8));
            
        end
    end
    
    
    
    % check convergence
    if mod(t,teval)==1
        
        conv = abs(mean(u(:))/mean(u_old(:))-1);
        
        if conv<tol
            break
        else
            u_old = u;
        end
    end
    
    
    for k=1:NPOP
        for j=1:NY
            for i=1:NX
                
                % Compute equilibrium distribution
                feq(i,j,k)=w(k)*rho(i,j)*(1 + 3*(u(i,j)*cx(k)+v(i,j)*cy(k))...
                    +9/2*(u(i,j)*cx(k)+v(i,j)*cy(k)).^2-3/2*(u(i,j)+v(i,j)).^2);
                
                Source(k)=w(k).*3.*rho(i,j)*(cx(k).*forcex+cy(k).*forcey);
                
                % Collision step (BGK collision operator)
                fcoll(i,j,k) = f(i,j,k) - 1/tau*(f(i,j,k)-feq(i,j,k))+Source(k);
            end
        end
    end
    
    % Streaming step (with periodic BC included)
    for j=1:NY
        
        if j>1
            jn = j-1;
        else
            jn = NY;
        end
        
        if j<NY
            jp=j+1;
        else
            jp=1;
        end
        
        
        for i=1:NX
            
            if i>1
                in=i-1;
            else
                in=NX;
            end
            
            if i<NX
                ip=i+1;
            else
                ip=1;
            end
            
            f(ip,j,1)  =  fcoll(i,j,1);
            f(i,jp,2)  =  fcoll(i,j,2);
            f(in,j,3)  =  fcoll(i,j,3);
            f(i,jn,4)  =  fcoll(i,j,4);
            f(ip,jp,5) =  fcoll(i,j,5);
            f(in,jp,6) =  fcoll(i,j,6);
            f(in,jn,7) =  fcoll(i,j,7);
            f(ip,jn,8) =  fcoll(i,j,8);
            f(i,j,9)   =  fcoll(i,j,9);
            
        end
    end
    
    
    
    % Boundary condition (NEBB)
    % Top wall (rest)
    rho(:,NY)=1/(1+v(:,NY))*(f(:,NY,9)+f(:,NY,1)+f(:,NY,3)+...
        2*(f(:,NY,2)+f(:,NY,5)+f(:,NY,6))-1/2.*forcey(:));
    u(:,NY)=0;
    v(:,NY)=0;
    
    f(:,NY,4)=f(:,NY,2)-2/3*v(:,NY)+1/6*forcey(:);
    f(:,NY,7)=f(:,NY,5)-1/6*v(:,NY)+0.5*(f(:,NY,1)-f(:,NY,3))-1/2*u(:,NY)+1/4*forcex(:)+1/6*forcey(:);
    f(:,NY,8)=f(:,NY,6)-1/6*v(:,NY)+0.5*(f(:,NY,3)-f(:,NY,1))+1/2*u(:,NY)-1/4*forcex(:)+1/6*forcey(:);
    
    % Bottom wall (rest)
    rho(:,1)=1/(1-v(:,1))*(f(:,1,9)+f(:,1,1)+f(:,1,3)+...
        2*(f(:,1,4)+f(:,1,7)+f(:,1,8))+1/2*forcey(:));
    u(:,1)=0;
    v(:,1)=0;
    
    f(:,1,2)=f(:,1,4)+2/3*v(:,1)-1/6*forcey(:);
    f(:,1,5)=f(:,1,7)+1/6*v(:,1)-0.5*(f(:,1,1)-f(:,1,3))+1/2*u(:,1)-1/4*forcex(:)-1/6*forcey(:);
    f(:,1,6)=f(:,1,8)+1/6*v(:,1)+0.5*(f(:,1,1)-f(:,1,3))-1/2*u(:,1)+1/4*forcex(:)-1/6*forcey(:);
    
    
    
end

% Half-force correction to recover physical fluid velocity
u=u+0.5*forcex;
v=v+0.5*forcex;


% Calculate performance information after the simulation is finished
runtime = toc(tstart);

% Compute error: L2 norm
error=zeros(NX,1);
for i=1:NX
    error(i)=(sqrt(sum((u(i,:)-u_analy).^2)))./sqrt(sum(u_analy.^2));
end
L2=1/NX*sum(error);

% Accuracy information
fprintf(' ----- accuracy information -----\n');
fprintf('        L2(u): %u\n',L2);
