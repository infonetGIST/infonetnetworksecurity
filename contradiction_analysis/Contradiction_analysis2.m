clc; clear;
%%
EbNodB = 5;

EbNo = 10.^(EbNodB./10);
p = 0.5*erfc(sqrt(2*EbNo)/2);

% p=0:0.01:0.5;
pa = 0:0.01:1;
% compromise rate
rc = 0.15;

% mode = 1; Unanimous rule
% mode = 2; Majority rule
mode = 2;

% [X, Y]=meshgrid(p,p_a);
% [W, T]=meshgrid(EbNodB,p_a);

[Z_comp, Z_usual] = regular_ldpc_analysis_fun3(p,pa,rc,mode);

figure(1);
plot(pa,Z_comp);

figure(2);
plot(pa,Z_usual);
% xlabel('p');ylabel('P(a)'); zlabel('Contradiction Probability');
% xlabel('EbNo');ylabel('P(a)'); zlabel('Contradiction Probability');
%%