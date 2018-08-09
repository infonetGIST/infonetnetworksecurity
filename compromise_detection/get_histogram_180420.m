clc;
% clear;
%%
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/LDPC_decoder');
%%
% Ns = 50; Nr = 100; dv = 10; dc = 6; num_avg = 186;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');
Ns = 100; Nr = 100; dv = 8; dc = 9; num_avg = 231;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
% Ns = 200; Nr = 100; dv = 6; dc = 13; num_avg = 317;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');
code_rate = Ns / (Ns+Nr);
%%
EsNodB = 5;
mode = 1; rc = 0.15; pa = 0.3;
num_trial = 1; 

for index=1:num_trial
    [relay_contradict] = wb_get_relay_contradict(Ns,Nr,rc,pa,EsNodB,num_avg);
    % load AveLRContradict_5dB AveLRContradict;
    %% K-means clustering
    ini_point = [0;num_avg];
    opts = statset('Display','off','MaxIter',100);
    [idx,C,sumd,D] = kmeans(relay_contradict,2,'Distance','sqeuclidean',...
    'Replicates',1,'Start',ini_point,'Options',opts);
    %     'Replicates',1,'Options',opts);    
    [a,b] = max(C);
    est_idx1 = find(idx==1);
    est_idx2 = find(idx==2);
    figure();
    plot(est_idx1,relay_contradict(est_idx1,1),'r.','MarkerSize',12); hold on;
    plot(est_idx2,relay_contradict(est_idx2,1),'b.','MarkerSize',12); hold off;
%%   
    usual_n = zeros(1,num_avg+1);
    for i=0:num_avg
        usual_n(i+1) = histc(relay_contradict(est_idx1,1),i);
    end
    figure();
    plot(0:1:num_avg,usual_n,'b-o');hold on;
%     axis([0 1 0 6]);

    comp_n = zeros(1,num_avg+1);
    for i=0:num_avg
        comp_n(i+1) = histc(relay_contradict(est_idx2,1),i);
    end
    plot(0:1:num_avg,comp_n,'r-v');
    
    xlabel('APC');ylabel('Number of Relays');    
end
%%
% rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix');
% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/LDPC_decoder');

rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');