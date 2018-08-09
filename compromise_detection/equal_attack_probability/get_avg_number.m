clc;clear;
%%
ini_num_avg = 112;
% num_avg = 110;
% pa = 0.1:0.1:1.0;
pa = 1.0;
% rc = 0.1:0.1:0.5;
rc = 0.25;

EbNodB = 5;

itr_num = 1;
num_node = 100;
%%
% D_Prob = zeros(length(pa),length(EbNodB));
% MD_Prob = zeros(length(pa),length(EbNodB));
% FA_Prob = zeros(length(pa),length(EbNodB));
Observe_Window = zeros(length(rc),length(pa));
for idx_rc = 1:length(rc)
    rrc = rc(idx_rc);
    fprintf('rc = %f\n',rrc);
    num_avg = ini_num_avg;
    for idx_pa = 1:length(pa)
        ppa = pa(idx_pa);
        fprintf('pa = %f\n',ppa);
        
        while 1
            fprintf('number of average = %d\n',num_avg);
            n_comp = ceil(num_node*rrc); n_usual = num_node - n_comp;            
            
            same = 0; MD = 0; FA = 0;
            for itr = 1:itr_num
                fprintf('iteration number = %d\n',itr);
                mode = 1;
                [AveLRContradict] = wb_Main_DAS_NetwrokCoding(rrc,ppa,EbNodB,num_avg,mode);                
        %     load AveLRContradict_5dB AveLRContradict;
                relay_contradict = AveLRContradict(101:200)';
        %% K-means clustering
                ini_point = [0;num_avg];
                opts = statset('Display','off','MaxIter',100);
                [idx,C,sumd,D] = kmeans(relay_contradict,2,'Distance','sqeuclidean',...
                'Replicates',1,'Start',ini_point,'Options',opts);
                %     'Replicates',1,'Options',opts);    
                [a,b] = max(C);
                est_idx1 = find(idx==1);
                est_idx2 = find(idx==2);
                figure;
                plot(est_idx1,relay_contradict(est_idx1,1),'r.','MarkerSize',12); hold on;
                plot(est_idx2,relay_contradict(est_idx2,1),'b.','MarkerSize',12); hold off;
        %%
                Knon_Attack_Posi = 1:1:n_comp;
                Knon_Attack_Posi = Knon_Attack_Posi';                
                AttackRelayIndex = find(idx==b);
                same = same + length(intersect(AttackRelayIndex,Knon_Attack_Posi));
                MD = MD + length(setdiff(Knon_Attack_Posi,AttackRelayIndex));
                FA = FA + length(setdiff(AttackRelayIndex,Knon_Attack_Posi));
            end
            D_Prob = same/(n_comp*itr);
            MD_Prob = MD/(n_comp*itr);
            FA_Prob = FA/((num_node - n_comp)*itr);
            
            fprintf('false alarm probability is = %f\n\n',FA_Prob);
            fprintf('miss detection probability is = %f\n\n',MD_Prob);
            if FA_Prob > 1e-3 || (1 - D_Prob) > 1e-3
                Observe_Window(idx_rc,idx_pa) = num_avg;
                num_avg = num_avg + 1;
                break;
            end
            num_avg = num_avg - 1;
        end
    end
end

% figure();
% plot(pa,Observe_Window(1,:), 'b'); hold on;
% plot(num_avg,MD_Prob, 'm');
% plot(num_avg,FA_Prob, 'r'); hold off;
% legend('Detection Probability','Miss Detection Probability','False Alarm Probability');
%%