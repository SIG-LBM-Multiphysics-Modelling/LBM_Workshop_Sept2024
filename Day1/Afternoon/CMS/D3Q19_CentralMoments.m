%--------------------------------------------------------------------------
% This program is part of the paper "?Multiphysics flow simulations at low
% cost using D3Q19 lattice Boltzmann methods based on central moments".
% Copyright (C) 2020  Alessandro De Rosis (derosis.alessandro@icloud.com)
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
% This is the "D3Q19_CentralMoments.m" file.
% Author: Alessandro De Rosis (derosis.alessandro@icloud.com)
% Day: 1st Jun 2020
%--------------------------------------------------------------------------

clear all
clc
syms U V W omega k5_star k6_star k7_star k8_star k9_star Fx Fy Fz R real

%% Define useful quantities
% Relaxation matrix
K = diag([1 1 1 1 1 omega omega omega omega omega 1 1 1 1 1 1 1 1 1]);
% Lattice directions and weights
w0 = 1/3;
ws = 1/18;
wl = 1/36;

cx = [0 1 -1 0 0 0 0 1 -1 1 -1 1 -1 1 -1 0 0 0 0];
cy = [0 0 0 1 -1 0 0 1 -1 -1 1 0 0 0 0 1 -1 1 -1];
cz = [0 0 0 0 0 1 -1 0 0 0 0 1 -1 -1 1 1 -1 -1 1];
w = [w0 ws ws ws ws ws ws wl wl wl wl wl wl wl wl wl wl wl wl];

T = sym(zeros(19,19));
M = zeros(19,19);
feq = sym(zeros(19,1));
U2 = U*U;
V2 = V*V;
W2 = W*W;
u2 = U*U+V*V+W*W;
for i=1:length(cx)
        feq(i) = w(i)*R*(1+3*(U*cx(i)+V*cy(i)+W*cz(i)) + 4.5*((U*cx(i)+V*cy(i)+W*cz(i))^2 - 1.5*(U*U+V*V+W*W)));
% 
%     % build equilibrium populations
%     if(cx(i)==0 && cy(i)==0 && cz(i)==0)
%         feq(i) = R/3*(1-(U2+V2+W2)+3*(U2*V2+U2*W2+V2*W2));
%     elseif(cx(i)~=0 && cy(i)==0 && cz(i)==0)
%         feq(i) = R/18*(1+3*cx(i)*U+3*(U2-V2-W2)-9*cx(i)*(U*V2+U*W2)-9*(U2*V2+U2*W2));
%     elseif(cx(i)==0 && cy(i)~=0 && cz(i)==0)
%         feq(i) = R/18*(1+3*cy(i)*V+3*(-U2+V2-W2)-9*cy(i)*(U2*V+V*W2)-9*(U2*V2+V2*W2));
%     elseif(cx(i)==0 && cy(i)==0 && cz(i)~=0)
%         feq(i) = R/18*(1+3*cz(i)*W+3*(-U2-V2+W2)-9*cz(i)*(U2*W+V2*W)-9*(U2*W2+V2*W2));
%     elseif(cx(i)~=0 && cy(i)~=0 && cz(i)==0)
%         feq(i) = R/36*(1+3*(cx(i)*U+cy(i)*V)+3*(U2+V2)+9*cx(i)*cy(i)*U*V+...
%                  9*(cy(i)*U2*V+cx(i)*U*V2)+9*U2*V2);
%     elseif(cx(i)~=0 && cy(i)==0 && cz(i)~=0)
%         feq(i) = R/36*(1+3*(cx(i)*U+cz(i)*W)+3*(U2+W2)+9*cx(i)*cz(i)*U*W+...
%                  9*(cz(i)*U2*W+cx(i)*U*W2)+9*U2*W2);
%     elseif(cx(i)==0 && cy(i)~=0 && cz(i)~=0)
%         feq(i) = R/36*(1+3*(cy(i)*V+cz(i)*W)+3*(V2+W2)+9*cy(i)*cz(i)*V*W+...
%                  9*(cz(i)*V2*W+cy(i)*V*W2)+9*V2*W2);
    %end
    
    % Set the basis
    CX = cx(i)-U;
    CY = cy(i)-V;
    CZ = cz(i)-W;
    CX2 = CX^2;
    CY2 = CY^2;
    CZ2 = CZ^2;
    T(1,i) = 1;
    T(2,i) = CX;
    T(3,i) = CY;
    T(4,i) = CZ;
    T(5,i) = CX2+CY2+CZ2;
    T(6,i) = CX2-CY2;
    T(7,i) = CY2-CZ2;
    T(8,i) = CX*CY;
    T(9,i) = CX*CZ;
    T(10,i) = CY*CZ;
    T(11,i) = CX2*CY;
    T(12,i) = CX*CY2;
    T(13,i) = CX2*CZ;
    T(14,i) = CX*CZ2;
    T(15,i) = CY2*CZ;
    T(16,i) = CY*CZ2;
    T(17,i) = CX2*CY2;
    T(18,i) = CX2*CZ2;
    T(19,i) = CY2*CZ2;
    
    CX = cx(i);
    CY = cy(i);
    CZ = cz(i);
    CX2 = CX^2;
    CY2 = CY^2;
    CZ2 = CZ^2;
    M(1,i) = 1;
    M(2,i) = CX;
    M(3,i) = CY;
    M(4,i) = CZ;
    M(5,i) = CX2+CY2+CZ2;
    M(6,i) = CX2-CY2;
    M(7,i) = CY2-CZ2;
    M(8,i) = CX*CY;
    M(9,i) = CX*CZ;
    M(10,i) = CY*CZ;
    M(11,i) = CX2*CY;
    M(12,i) = CX*CY2;
    M(13,i) = CX2*CZ;
    M(14,i) = CX*CZ2;
    M(15,i) = CY2*CZ;
    M(16,i) = CY*CZ2;
    M(17,i) = CX2*CY2;
    M(18,i) = CX2*CZ2;
    M(19,i) = CY2*CZ2;

 end
T = simplify(T);
N = simplify(T*M^(-1)); %shift matrix
% N = [1,        0,        0,        0,             0,                 0,                 0,     0,     0,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     -U,        1,        0,        0,             0,                 0,                 0,     0,     0,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     -V,        0,        1,        0,             0,                 0,                 0,     0,     0,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     -W,        0,        0,        1,             0,                 0,                 0,     0,     0,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     U^2 + V^2 + W^2,     -2*U,     -2*V,     -2*W,             1,                 0,                 0,     0,     0,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     U^2 - V^2,     -2*U,      2*V,        0,             0,                 1,                 0,     0,     0,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     V^2 - W^2,        0,     -2*V,      2*W,             0,                 0,                 1,     0,     0,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     U*V,       -V,       -U,        0,             0,                 0,                 0,     1,     0,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     U*W,       -W,        0,       -U,             0,                 0,                 0,     0,     1,     0,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     V*W,        0,       -W,       -V,             0,                 0,                 0,     0,     0,     1,    0,    0,    0,    0,    0,    0, 0, 0, 0;...
%     -U^2*V,    2*U*V,      U^2,        0,          -V/3,          -(2*V)/3,              -V/3,  -2*U,     0,     0,    1,    0,    0,    0,    0,    0, 0, 0, 0;...
%     -U*V^2,      V^2,    2*U*V,        0,          -U/3,               U/3,              -U/3,  -2*V,     0,     0,    0,    1,    0,    0,    0,    0, 0, 0, 0;...
%     -U^2*W,    2*U*W,        0,      U^2,          -W/3,          -(2*W)/3,              -W/3,     0,  -2*U,     0,    0,    0,    1,    0,    0,    0, 0, 0, 0;...
%     -U*W^2,      W^2,        0,    2*U*W,          -U/3,               U/3,           (2*U)/3,     0,  -2*W,     0,    0,    0,    0,    1,    0,    0, 0, 0, 0;...
%     -V^2*W,        0,    2*V*W,      V^2,          -W/3,               W/3,              -W/3,     0,     0,  -2*V,    0,    0,    0,    0,    1,    0, 0, 0, 0;...
%     -V*W^2,        0,      W^2,    2*V*W,          -V/3,               V/3,           (2*V)/3,     0,     0,  -2*W,    0,    0,    0,    0,    0,    1, 0, 0, 0;...
%     U^2*V^2, -2*U*V^2, -2*U^2*V,        0, U^2/3 + V^2/3, (2*V^2)/3 - U^2/3,     U^2/3 + V^2/3, 4*U*V,     0,     0, -2*V, -2*U,    0,    0,    0,    0, 1, 0, 0;...
%     U^2*W^2, -2*U*W^2,        0, -2*U^2*W, U^2/3 + W^2/3, (2*W^2)/3 - U^2/3, W^2/3 - (2*U^2)/3,     0, 4*U*W,     0,    0,    0, -2*W, -2*U,    0,    0, 0, 1, 0;...
%     V^2*W^2,        0, -2*V*W^2, -2*V^2*W, V^2/3 + W^2/3,   - V^2/3 - W^2/3, W^2/3 - (2*V^2)/3,     0,     0, 4*V*W,    0,    0,    0,    0, -2*W, -2*V, 0, 0, 1];
% N = simplify(N);
%N = eye(19,19);
%% Compute central moments
T = simplify(T);
%T = M;
K_eq = simplify(T*feq);
Id = eye(19,19);
K_pre = sym(zeros(19,1));
syms k5_pre k6_pre k7_pre k8_pre k9_pre k24_pre k25_pre k26_pre real
K_pre(1) = R;
K_pre(6) = k5_pre;
K_pre(7) = k6_pre;
K_pre(8) = k7_pre;
K_pre(9) = k8_pre;
K_pre(10) = k9_pre;
%post-collision central moments
K_star = (Id-K)*K_pre + K*K_eq
%post collision populations
syms k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 real
K_sym = [R k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18];
for i=1:19
    if(K_star(i)~=sym(0))
        K_star(i) = K_sym(i);
    end
end
raw_moments = collect(simplify(N^(-1)*K_star), K_star)
syms r0 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 r14 r15 r16 r17 r18 real
r = [r0 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 r14 r15 r16 r17 r18]'; %symbolic raw moments
f_post_collision_twosteps = simplify(M\r)