clc;clear;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis\');
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/LDPC_decoder');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix');
%%
num_avg = 10:10:350;
% num_avg = 447;
pa = 0.1; rc = 0.5;
% pa = 1; rc = 0.1;
% num_node = 100;

% Ns = 50; Nr = 100; dv = 10; dc = 6;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');
% Ns = 100; Nr = 100; dv = 8; dc = 9;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
Ns = 200; Nr = 100; dv = 6; dc = 13;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');

code_rate = Ns / (Ns+Nr);
% EbNodB = 8;
% EbNo = 10.^(EbNodB./10);
% EsNo = EbNo * code_rate;
% EsNodB = 10*log10(EsNo);
EsNodB = 5;
EsNo = 10.^(EsNodB./10);

itr_num = 20;
% itr_num = 1;
n_comp = ceil(Nr*rc); n_usual = Nr - n_comp;
%%
D_Prob = zeros(length(num_avg),1);
MD_Prob = zeros(length(num_avg),1);
FA_Prob = zeros(length(num_avg),1);
tic
for ind = 1:length(num_avg)
    nnum_avg = num_avg(ind);
    fprintf('Observation window = %d\n',nnum_avg);
    same = 0; MD = 0; FA = 0;
    for itr = 1:itr_num
        fprintf('iteration number = %d\n',itr);
        mode = 1;            
%                 [AveLRContradict] = wb_get_AveLRContradict(rrc,ppa,EbNodB,num_avg,mode);
        [AveLRContradict,tmp_pos] = wb_Main_DAS_NetwrokCoding(Ns,Nr,rc,pa,EsNodB,nnum_avg,mode);
        %     load AveLRContradict_5dB AveLRContradict;
        relay_contradict = AveLRContradict(Ns+1:end)';
        %% K-means clustering
        ini_point = [0;nnum_avg];
        opts = statset('Display','off','MaxIter',100);
        [idx,C,sumd,D] = kmeans(relay_contradict,2,'Distance','sqeuclidean',...
           'Replicates',1,'Start',ini_point,'Options',opts);
       [a,b] = max(C);
       est_idx1 = find(idx==1);
       est_idx2 = find(idx==2);
%                 figure;
%                 plot(est_idx1,relay_contradict(est_idx1,1),'r.','MarkerSize',12); hold on;
%                 plot(est_idx2,relay_contradict(est_idx2,1),'b.','MarkerSize',12); hold off;
        %%
%        Knon_Attack_Posi = 1:1:n_comp;
       Knon_Attack_Posi = tmp_pos;
       Knon_Attack_Posi = Knon_Attack_Posi';                
       AttackRelayIndex = find(idx==b);
       same = same + length(intersect(AttackRelayIndex,Knon_Attack_Posi));
       MD = MD + length(setdiff(Knon_Attack_Posi,AttackRelayIndex));
       FA = FA + length(setdiff(AttackRelayIndex,Knon_Attack_Posi));
    end
    D_Prob(ind,1) = same/(n_comp*itr)*100;
    MD_Prob(ind,1) = MD/(n_comp*itr)*100;
    FA_Prob(ind,1) = FA/((Nr - n_comp)*itr)*100;
end
toc
figure();
plot(num_avg,D_Prob, 'b'); hold on;
% plot(num_avg,MD_Prob, 'm');
plot(num_avg,FA_Prob, 'r'); hold off;
% legend('Detection Probability','Miss Detection Probability','False Alarm Probability');
legend('Detection Probability','False Alarm Probability');
grid on;
ax = gca;
ax.GridAlpha = 0.5;
ax.GridColor = [0.1 0.1 0.1];
ax.MinorGridAlpha = 0.5;
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis\');
% rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/LDPC_decoder');
% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix');

% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');

rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');
