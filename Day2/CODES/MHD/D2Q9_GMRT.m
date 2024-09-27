%--------------------------------------------------------------------------
% This program is part of the paper "﻿Double-D2Q9 lattice Boltzmann models
% with extended equilibrium for two-dimensional magnetohydrodynamics
% flows".
% Copyright (C) 2020  Alessandro De Rosis
% (alessandro.derosis@manchester.ac.uk,
% derosis.alessandro@icloud.com)
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% This is the "﻿D2Q9_GMRT_MHD.m" file.
% Author: Alessandro De Rosis (alessandro.derosis@manchester.ac.uk,
% derosis.alessandro@icloud.com)
% Day: 1st Dec 2020
%--------------------------------------------------------------------------
clear all
clc
addpath('/Applications/Matlab_R2019b.app/toolbox/symbolic/symbolic')
%% Initialize some symbolic variables
syms U V R omega f0 f1 f2 f3 f4 f5 f6 f7 f8 BX BY...
     k4_pre k5_pre...
     k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star...
     r0 r1 r2 r3 r4 r5 r6 r7 r8 real
%% Define lattice directions, weights and other useful quantities of the D2Q9 model
raw_moments = false; % if false -> central moments
cx = [0 1 0 -1 0 1 -1 -1 1];
cy = [0 0 1 0 -1 1 1 -1 -1];
w = [4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.];
np = length(cx);
cs = 1/sqrt(3);
cs2 = cs^2;
cs3 = cs^3;
cs4 = cs^4;
cs6 = cs^6;
cs8 = cs^8;
f = [f0 f1 f2 f3 f4 f5 f6 f7 f8]';
feq = sym(zeros(np,1));
% transformation matrices
T = sym(zeros(np,np)); 
M = zeros(np,np); 
% relaxation matrix
Lambda = diag([1, 1, 1, 1, omega, omega, 1, 1, 1]);
% identity matrix
Id = eye(np,np); 
for i=1:np
    % build the complete equilibria
    first_order = 1/cs2*(U*cx(i)+V*cy(i));
    second_order = 1/(2*cs4)*((cx(i)*cx(i)-1/3)*U^2+...
                              (cy(i)*cy(i)-1/3)*V^2+...
                              2*cx(i)*cy(i)*U*V);
    third_order = 1/(2*cs6)*((cx(i)^2-1/3)*cy(i)*U*U*V+(cy(i)^2-1/3)*cx(i)*U*V*V);
    fourth_order = 1/(4*cs8)*((cx(i)^2-1/3)*(cy(i)^2-1/3)*U*U*V*V);
    termMHD = w(i)/(2*cs4)*(0.5*(BX*BX+BY*BY)*(cx(i)*cx(i)+cy(i)*cy(i))-(cx(i)*BX+cy(i)*BY)^2);
    feq(i) = w(i)*R*(1+first_order+second_order+third_order+fourth_order)+termMHD;
    
    % build the transformation matrix T 
    CX = cx(i)-U;
    CY = cy(i)-V;
    T(1,i) = 1;
    T(2,i) = CX;
    T(3,i) = CY;
    T(4,i) = CX*CX+CY*CY;
    T(5,i) = CX*CX-CY*CY;
    T(6,i) = CX*CY;
    T(7,i) = CX*CX*CY;
    T(8,i) = CX*CY*CY;
    T(9,i) = CX*CX*CY*CY;
    
    % build the tranformation matrix M
    CX = cx(i);
    CY = cy(i);
    M(1,i) = 1;
    M(2,i) = CX;
    M(3,i) = CY;
    M(4,i) = CX*CX+CY*CY;
    M(5,i) = CX*CX-CY*CY;
    M(6,i) = CX*CY;
    M(7,i) = CX*CX*CY;
    M(8,i) = CX*CY*CY;
    M(9,i) = CX*CX*CY*CY;
end
T = simplify(T);
if(raw_moments==true)
    T = M;
end
% shift matrix
N = simplify(T*M^(-1)); 
%% Collision stage
%pre-collision central moments
K_pre = simplify(T*f); 
% equilibrium central moments
K_eq = simplify(T*feq); 
K_pre(5) = k4_pre;
K_pre(6) = k5_pre;
% post-collision central moments
K_star = simplify((Id-Lambda)*K_pre + Lambda*K_eq) 
% post-collision populations
K_sym = [R k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star];
for i=1:np
    if(K_star(i)~=sym(0))
        K_star(i) = K_sym(i);
    end
end
raw_moments = simplify(N^(-1)*K_star)
r = [r0 r1 r2 r3 r4 r5 r6 r7 r8]';
f_post_collision = collect(simplify(M^(-1)*r),K_star)