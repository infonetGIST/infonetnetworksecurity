% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: detection_eval.m
% this script generates detection and false alarm probabilities with varying observation window
% *****************************************************************
%% Parameters
addpath('../H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');
addpath('../LDPC_decoder');
% case 1
Ns = 50; Nr = 100; dv = 10; dc = 6; 
addpath('../H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');
% case 2
% Ns = 100; Nr = 100; dv = 8; dc = 9;
% addpath('../\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
% case 3
% Ns = 200; Nr = 100; dv = 6; dc = 13;
% addpath('../\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');

code_rate = Ns / (Ns+Nr);
EsNodB = 5;
mode = 1; rc = 0.15; pa = 0.3;
num_trial = 1; 
num_avg = 203;
for index=1:num_trial
    [relay_contradict] = wb_Main_DAS_NetwrokCoding(Ns,Nr,rc,pa,EsNodB,num_avg,mode);
    relay_contradict = relay_contradict';
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
%% Draw histogram
    usual_n = zeros(1,num_avg+1);
    for i=0:num_avg
        usual_n(i+1) = histc(relay_contradict(est_idx1,1),i);
    end
    figure(); plot(0:1:num_avg,usual_n,'b-o');hold on;

    comp_n = zeros(1,num_avg+1);
    for i=0:num_avg
        comp_n(i+1) = histc(relay_contradict(est_idx2,1),i);
    end
    plot(0:1:num_avg,comp_n,'r-v');
    
    xlabel('APC');ylabel('Number of Relays');    
end
%%
rmpath(genpath('../'));