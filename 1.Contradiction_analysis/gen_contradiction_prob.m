% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: gen_contradiction_prob.m
% this script generates contradiction probabilities of compromised relays and usual relays
% *****************************************************************
%% Parameters
% Ns: # of sources, Nr: # of relays, dv: # of edges per variable node, dc: # of edges per check node
% case 1
% Ns = 50; Nr = 100; dv = 10; dc = 6; 
% case 2
Ns = 100; Nr = 100; dv = 8; dc = 9;
% case 3
% Ns = 200; Nr = 100; dv = 6; dc = 13; 
EsNodB = -4:.1:10;
pa = 0:0.01:1; % attack probability
rc = 0.05; % compromise rate
mode = 1; % mode = 1; Unanimity rule, mode = 2; Majority rule
%%
addpath('../MP_analysis');
EsNo = 10.^(EsNodB./10);

p = 0.5*erfc(sqrt(EsNo)); % channel error probability

[X, Y] = meshgrid(p,pa);
[W, T] = meshgrid(EsNodB,pa);

p1 = X+Y-2.*X.*Y; % attack probability of the compromised relay

% iterative message passing decoding
% pe_m : probability that even number of errors in (dc - 1) source nodes
% pe_m : probability that odd number of errors in (dc - 1) source nodes
[pe_m, po_m, p_out] = regular_ldpc_analysis_fun(p,p1,rc,dv,dc,mode);

% For compromised node,
Z = pe_m.*(p1) + po_m.*(1-(p1)); % contradiction probability of compromised relay node
% For usual node,
Z2 = pe_m.*(X) + po_m.*(1-(X)); % contradiction probability of usual relay node
%% Draw figures
figure(1)
mesh(W,Y,Z)
xlabel('E_s/N_o [dB]');ylabel('P_a'); zlabel('Contradiction Probability');
figure(2)
mesh(W,Y,Z2)
xlabel('E_s/N_o [dB]');ylabel('P_a'); zlabel('Contradiction Probability');
%%
rmpath(genpath('../'));