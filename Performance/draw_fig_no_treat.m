clc;clear;
%%
addpath(genpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\Performance\No_treatment'));
% addpath(genpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/Performance/No_treatment'));
%%
load BERRecord_nocoding_awgn.mat;
load BERRecord_noattack_w_coding.mat
load BERRecord_rc05_pa50.mat;
load BERRecord_rc05_pa100.mat;
load BERRecord_rc15_pa50.mat;
load BERRecord_rc15_pa100.mat;
%%
figure;
semilogy(0:1:10, BERRecord_nocoding_awgn); hold on;
semilogy(0:1:5, BERRecord_noattack_w_coding);
semilogy(0:1:9, BERRecord_rc05_pa50);
semilogy(0:1:9, BERRecord_rc05_pa100);
semilogy(0:1:9, BERRecord_rc15_pa50); 
semilogy(0:1:9, BERRecord_rc15_pa100);
hold off;
grid on; 
legend('No coding','No attack','Under attack R_c=0.05, p_a=0.5','Under attack R_c=0.05, p_a=1.0',...
    'Under attack R_c=0.15, p_a=0.5', 'Under attack R_c=0.15, p_a=1.0');
axis([0 9 1e-5 1]);
%%
rmpath(genpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\Performance\No_treatment'));
% rmpath(genpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/Performance/No_treatment'));