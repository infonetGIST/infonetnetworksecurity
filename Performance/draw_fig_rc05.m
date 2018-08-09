clc;clear;
%%
addpath(genpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\Performance\No_treatment'));
addpath(genpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\Performance\CD_AC'));
% addpath(genpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/Performance/No_treatment'));
% addpath(genpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/Performance/CD_AC'));
%%
load BERRecord_nocoding_awgn.mat;
load BERRecord_noattack_w_coding.mat
load BERRecord_rc05_pa100.mat;
load BERRecord_discard_rc05_pa50.mat;
load BERRecord_discard_rc05_pa100.mat;
load BERRecord_reverse_rc05_pa50.mat;
load BERRecord_reverse_rc05_pa100.mat;
%%
figure;
semilogy(0:1:10, BERRecord_nocoding_awgn); hold on;
semilogy(0:1:5, BERRecord_noattack_w_coding);
semilogy(0:1:9, BERRecord_rc05_pa100);
semilogy(0:1:5, BERRecord_discard_rc05_pa50);
semilogy(0:1:6, BERRecord_discard_rc05_pa100);
semilogy(0:1:9, BERRecord_reverse_rc05_pa50);
semilogy(0:1:5, BERRecord_reverse_rc05_pa100); 
hold off;
legend('No coding','No attack','Under attack','Discard, R_c=0.05, p_a=0.5',...
    'Discard, R_c=0.05, p_a=1.0','Reverse, R_c=0.05, p_a=0.5', 'Reverse, R_c=0.05, p_a=1.0');
axis([0 9 5e-6 1e0]);
grid on;
ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor = [0.1 0.1 0.1];
ax.MinorGridAlpha = 0.5;
%%
rmpath(genpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\Performance\No_treatment'));
rmpath(genpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\Performance\CD_AC'));
% rmpath(genpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/Performance/No_treatment'));
% rmpath(genpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/Performance/CD_AC'));