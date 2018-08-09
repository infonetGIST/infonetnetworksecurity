% -------------------------------------------------------------------------
% filename : draw_hisotram.m
% Objectives : this script generates detection and false alarm probabilities with varying observation window
% written by Woong-Bi Lee, July 2018
% -------------------------------------------------------------------------
%%
addpath('../MP_analysis');

load histo_180420.mat;

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
rmpath(genpath('../'));