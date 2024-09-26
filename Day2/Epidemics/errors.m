clc

A = load('dataBGK.txt');
BGK.t = A(:,1);
BGK.N = A(:,2);
BGK.S = A(:,3);
BGK.E = A(:,4);
BGK.I = A(:,5);
BGK.R = A(:,6);
BGK.D = A(:,7);

C = load('dataFD.txt');
FD.t = C(:,1);
FD.N = C(:,2);
FD.S = C(:,3);
FD.E = C(:,4);
FD.I = C(:,5);
FD.R = C(:,6);
FD.D = C(:,7);

error_L2_BGK_FD.N = norm(FD.N-BGK.N)/norm(FD.N)*100;
error_L2_BGK_FD.S = norm(FD.S-BGK.S)/norm(FD.S)*100;
error_L2_BGK_FD.E = norm(FD.E-BGK.E)/norm(FD.E)*100;
error_L2_BGK_FD.I = norm(FD.I-BGK.I)/norm(FD.I)*100;
error_L2_BGK_FD.R = norm(FD.R-BGK.R)/norm(FD.R)*100;
error_L2_BGK_FD.D = norm(FD.D-BGK.D)/norm(FD.D)*100;

error_Linf_BGK_FD.N = norm(FD.N-BGK.N, inf)/norm(FD.N, inf)*100;
error_Linf_BGK_FD.S = norm(FD.S-BGK.S, inf)/norm(FD.S, inf)*100;
error_Linf_BGK_FD.E = norm(FD.E-BGK.E, inf)/norm(FD.E, inf)*100;
error_Linf_BGK_FD.I = norm(FD.I-BGK.I, inf)/norm(FD.I, inf)*100;
error_Linf_BGK_FD.R = norm(FD.R-BGK.R, inf)/norm(FD.R, inf)*100;
error_Linf_BGK_FD.D = norm(FD.D-BGK.D, inf)/norm(FD.D, inf)*100;
