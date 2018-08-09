clc; clear;
%%
% Ns = 50; Nr = 100; dv = 10; dc = 6;
Ns = 100; Nr = 100; dv = 8; dc = 9;
% Ns = 200; Nr = 100; dv = 6; dc = 13;

code_rate = Ns / (Ns+Nr);
% EbNodB=[0 3 5];

EbNodB = -4:.1:10;
EbNo = 10.^(EbNodB./10);
EsNo = EbNo * code_rate;
% EsNodB = 0:.1:10;
% EsNo = 10.^(EsNodB./10);

p = 0.5*erfc(sqrt(EsNo));

% p=0:0.01:0.5;
p_a = 0:0.01:1;
rc = 0.05; % compromise rate

% mode = 1; Unanimous rule, mode = 2; Majority rule
mode = 1;

[X, Y]=meshgrid(p,p_a);
[W, T]=meshgrid(EbNodB,p_a);

p1 = X+Y-2.*X.*Y;

[pe_m, po_m, p_out]=regular_ldpc_analysis_fun4(p,p1,rc,dv,dc,mode);

% For attacked node,
Z = pe_m.*(p1) + po_m.*(1-(p1));
% For non-attacked node,
Z2 = pe_m.*(X) + po_m.*(1-(X));

% Z = p_out;
figure()
mesh(W,Y,Z)
% xlabel('p');ylabel('P(a)'); zlabel('Contradiction Probability');
xlabel('E_b/N_o [dB]');ylabel('P_a'); zlabel('Contradiction Probability');
figure();
mesh(W,Y,Z2)
% xlabel('p');ylabel('P(a)'); zlabel('Contradiction Probability');
xlabel('E_b/N_o [dB]');ylabel('P_a'); zlabel('Contradiction Probability');
%%