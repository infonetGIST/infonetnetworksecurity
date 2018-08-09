clc;clear;
% Rc = 0.05
% num_avg = 840;
% Rc = 0.15
% num_avg = 906;
num_avg = 1000;
EbNodB = 5;
mode = 1;
pc = 0.5; pa = 1.0;
num_trial = 1;

for index=1:num_trial
    [AveLRContradict] = wb_Main_DAS_NetwrokCoding(pc,pa,EbNodB,num_avg,mode);
    % load AveLRContradict_5dB AveLRContradict;
    relay_contradict = AveLRContradict(101:end)';
    EbNo=10.^(EbNodB./10);
    p=0.5*erfc(sqrt(2*EbNo)/2);
    % pa = 1.0; pc = 1.0;
    % [Prob_comp Prob_usual]=regular_ldpc_analysis_fun2(p,pa,pc);
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
%    %%
%     Knon_Attack_Posi = 1:1:n_comp;
%     Knon_Attack_Posi = Knon_Attack_Posi';                
%     AttackRelayIndex = find(idx==b);
%     same = same + length(intersect(AttackRelayIndex,Knon_Attack_Posi));
%     MD = MD + length(setdiff(Knon_Attack_Posi,AttackRelayIndex));
%     FA = FA + length(setdiff(AttackRelayIndex,Knon_Attack_Posi));
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
    
%     plot(Optimal_T,0:0.001:8,'g--'); 
    xlabel('APC');ylabel('Number of Relays');
%     legend('APC Histogram','Threshold');
    

% 
%     index_large_real = find(relay_contradict/num_avg > Optimal_T);   
    
    
end
hold off;