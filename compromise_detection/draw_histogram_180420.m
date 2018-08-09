clc;clear;
%%
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis');
%%
load histo_180420.mat;

N = 231;
EsNodB = 5; 
EsNo = 10.^(EsNodB./10);
p = 0.5*erfc(sqrt(EsNo));
dv = 8; dc = 9;
theo_pa = 0.0; theo_rc = 0.0;
mode_a = 1; %  mode = 1; Unanimous rule, mode = 2; Majority rule
[min_Prob_comp, min_Prob_usual]=regular_ldpc_analysis_fun3(p,theo_pa,theo_rc,dv,dc,mode_a);
theo_pa = 1.0; theo_rc = 0.0000000001;
[max_Prob_comp, max_Prob_usual]=regular_ldpc_analysis_fun3(p,theo_pa,theo_rc,dv,dc,mode_a);
min_APC_usual = N*min_Prob_usual;
max_APC_usual = N*max_Prob_usual;
min_APC_comp = N*min_Prob_comp;
max_APC_comp = N*max_Prob_comp;

figure();
subplot(3,1,1);
plot(0:1:num_avg,usual_n_rc15_pa30,'b-o');hold on;
plot(0:1:num_avg,comp_n_rc15_pa30,'r-v');
xlabel('APC');ylabel('Number of Relays');
% stem(min_APC_usual,7,'m');
% stem(max_APC_usual,7,'c');
% stem(max_APC_comp,7,'c');
axis([0 N 0 20]);

subplot(3,1,2);
plot(0:1:num_avg,usual_n_rc15_pa50,'b-o');hold on;
plot(0:1:num_avg,comp_n_rc15_pa50,'r-v');
xlabel('APC');ylabel('Number of Relays');
% stem(min_APC_usual,7,'m');
% stem(max_APC_usual,7,'c');
% stem(max_APC_comp,7,'c');
axis([0 N 0 10]);

subplot(3,1,3);
plot(0:1:num_avg,usual_n_rc15_pa100,'b-o');hold on;
plot(0:1:num_avg,comp_n_rc15_pa100,'r-v');
xlabel('APC');ylabel('Number of Relays');
% stem(min_APC_usual,7,'m');
% stem(max_APC_usual,7,'c');
% stem(max_APC_comp,7,'c');
axis([0 N 0 10]);
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis');