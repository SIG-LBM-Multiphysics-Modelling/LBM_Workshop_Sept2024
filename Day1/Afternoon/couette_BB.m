% This code is designed for the LBM Workshop 2024 @ Universtiy of
% Manchester. 
% This code is provided under the MIT license. 
% Author: Goncalo Silva
% 
% =========================================================================
% Example matlab code for computing a Couette flow with BB boundary


clear all
close all
clc

% simulation parameters
NX=3;             % channel length
NY=32;             % channel width
NSTEPS=1e6;     % number of simulation time steps



rho0=1;                 % rest density state
u_max=0.1;              % maximum velocity
Re=0.1;                 % Reynolds number; scaling parameter in simulation (Re=NY*u_max/nu)
nu=NY*u_max/Re;         % kinematic shear viscosity
tau=3*nu+1/2;           % relaxation time (BGK model)


% Lattice parameters; note zero direction is last
NPOP=9;                                         % number of velocities
w  = [1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9]; % weights
cx = [1 0 -1  0 1 -1 -1  1 0];                  % velocities, x components
cy = [0 1  0 -1 1  1 -1 -1 0];                  % velocities, y components


% Node locations
x = (1:NX)-0.5;
y = (1:NY)-0.5;

% Analytical solution: Couette velocity
u_analy=u_max/NY.*y;


% initialize populations
feq=zeros(NX,NY,NPOP);
for k=1:NPOP
        feq(:,:,k)=w(k); % assuming density equal one and zero velocity initial state
end

f=feq;
fcoll=feq;

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
                
                % Collision step (BGK collision operator)
                fcoll(i,j,k) = f(i,j,k) - 1/tau*(f(i,j,k)-feq(i,j,k));
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


    
    % Boundary condition (bounce-back)
    % Top wall (moving with tangential velocity u_max)
    f(:,NY,4)=fcoll(:,NY,2);
    f(:,NY,7)=fcoll(:,NY,5)-(1/6)*u_max;
    f(:,NY,8)=fcoll(:,NY,6)+(1/6)*u_max;
    
    % Bottom wall (rest)
    f(:,1,2)=fcoll(:,1,4);
    f(:,1,5)=fcoll(:,1,7);
    f(:,1,6)=fcoll(:,1,8);
    
    
    
end


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
