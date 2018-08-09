clc; clear;
%%
EbNodB = 5;

EbNo = 10.^(EbNodB./10);
p = 0.5*erfc(sqrt(2*EbNo)/2);

% p=0:0.01:0.5;
pa = 1.0;
% compromise rate
rc_dic = 0:0.01:1;

% mode = 1; Unanimous rule
% mode = 2; Majority rule
mode = 1;
dv = 4; dc = 5;
% [X, Y]=meshgrid(p,p_a);
% [W, T]=meshgrid(EbNodB,p_a);
for i = 1: length(rc_dic)
    [Z_comp(i), Z_usual(i)] = regular_ldpc_analysis_fun4_180411(p,pa,rc_dic(i),dv,dc,mode);
end
figure(1);
plot(rc_dic,Z_comp);

figure(2);
plot(rc_dic,Z_usual);
% xlabel('p');ylabel('P(a)'); zlabel('Contradiction Probability');
% xlabel('EbNo');ylabel('P(a)'); zlabel('Contradiction Probability');
%%