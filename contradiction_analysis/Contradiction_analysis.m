clc; clear;
%%
Ns = 50; Nr = 100;
code_rate = Ns / (Ns+Nr);
% EbNodB=[0 3 5];
% EbNodB = 0:.1:10;
% EbNo = 10.^(EbNodB./10);
% EsNo = EbNo * code_rate;
EsNodB = 0:.1:10;
EsNo = 10.^(EsNodB./10);

p = 0.5*erfc(sqrt(EsNo));

% p=0:0.01:0.5;
p_a = 0:0.01:1;
% compromise rate
rc = 0.5;

% mode = 1; Unanimous rule, mode = 2; Majority rule
mode = 1;
dv = 10; dc = 6;
% dv = 8; dc = 9;
% dv = 6; dc = 13;

[X, Y]=meshgrid(p,p_a);
[W, T]=meshgrid(EsNodB,p_a);

p1 = X+Y-2.*X.*Y;
sum_pe_m=0;
sum_po_m=0;

for i=1:1:1
    [pe_m, po_m, p_out]=regular_ldpc_analysis_fun2(X,p1,rc,dv,dc,mode);

    sum_pe_m = pe_m + sum_pe_m;
    sum_po_m = po_m + sum_po_m;
end

av_pe_m=sum_pe_m/i;
av_po_m=sum_po_m/i;

% error, even number of errors -> flip
% no error, odd number of errors -> flip
% Flipping probability

% {error occured, but flipped} or {error not occured, but flipped}
% For attacked node,
Z = av_pe_m.*(p1) + av_po_m.*(1-(p1));
% For non-attacked node,
Z2 = av_pe_m.*(X) + av_po_m.*(1-(X));

% Z = p_out;
figure(1)
mesh(W,Y,Z)
% xlabel('p');ylabel('P(a)'); zlabel('Contradiction Probability');
xlabel('E_b/N_o [dB]');ylabel('P_a'); zlabel('Contradiction Probability');
figure(2);
mesh(W,Y,Z2)
% xlabel('p');ylabel('P(a)'); zlabel('Contradiction Probability');
xlabel('E_b/N_o [dB]');ylabel('P_a'); zlabel('Contradiction Probability');
%%