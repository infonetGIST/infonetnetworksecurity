clc;clear;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis\');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
%%
EbNodB = 0:5;
% EbNodB = [0,5];
num_avg = 100;
itr_max = 10;
%%
D_Prob = zeros(length(EbNodB),1);
MD_Prob = zeros(length(EbNodB),1);
FA_Prob = zeros(length(EbNodB),1);
for IEsNO2 = 1:length(EbNodB)
    same = 0; MD = 0; FA = 0;
    for itr = 1:itr_max     
        EbNo = 10.^(EbNodB(IEsNO2)./10);
        p = 0.5*erfc(sqrt(2*EbNo)/2);
        theo_pa = 0.0; theo_rc = 0.0;
        mode_a = 2; %  mode = 1; Unanimous rule, mode = 2; Majority rule
        [Prob_comp min_Prob_usual]=regular_ldpc_analysis_fun3(p,theo_pa,theo_rc,mode_a);
        theo_pa = 1.0; theo_rc = 1.0;
        [Prob_comp max_Prob_usual]=regular_ldpc_analysis_fun3(p,theo_pa,theo_rc,mode_a);

        % usual relays
        p2 = max_Prob_usual; % contradiction probability of usual relays
        m2 = num_avg*p2; s2 = sqrt(num_avg*p2*(1-p2));
        y2 = pdf('normal',0:num_avg,m2,s2);

        mode = 1;
        rc = 0.15; pa = 1.0;
        n_comp = ceil(num_avg*rc);
        n_usual = num_avg - n_comp;
        [AveLRContradict] = wb_get_AveLRContradict(rc,pa,EbNodB(IEsNO2),num_avg,mode);
    %     load AveLRContradict_5dB AveLRContradict;

        relay_contradict = AveLRContradict(101:200)';
    %% K-means clustering
        ini_point = [0;100];
        opts = statset('Display','final','MaxIter',100);
        [idx,C,sumd,D] = kmeans(relay_contradict,2,'Distance','cityblock',...
        'Replicates',1,'Start',ini_point,'Options',opts);
        %     'Replicates',1,'Options',opts);    
        est_idx1 = find(idx==1);
        est_idx2 = find(idx==2);

    %     figure;
    %     plot(est_idx1,relay_contradict(est_idx1,1),'r.','MarkerSize',12); hold on;
    %     plot(est_idx2,relay_contradict(est_idx2,1),'b.','MarkerSize',12); hold off;
    %%
        Knon_Attack_Posi = 1:1:n_comp;
        AttackRelayIndex = est_idx2;
        same = same + length(intersect(AttackRelayIndex,Knon_Attack_Posi));
        MD = MD + length(setdiff(Knon_Attack_Posi,AttackRelayIndex));
        FA = FA + length(setdiff(AttackRelayIndex,Knon_Attack_Posi));
    end
    D_Prob(IEsNO2) = same/(n_comp*itr_max);
    MD_Prob(IEsNO2) = MD/(n_comp*itr_max);
    FA_Prob(IEsNO2) = FA/(n_usual*itr_max);    
end

figure();
plot(EbNodB,D_Prob, 'b'); hold on;
plot(EbNodB,MD_Prob, 'm');
plot(EbNodB,FA_Prob, 'r'); hold off;
legend('Detection Probability','Miss Detection Probability','False Alarm Probability');

%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis\');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');