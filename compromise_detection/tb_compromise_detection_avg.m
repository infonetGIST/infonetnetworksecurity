clc;clear;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis\');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
%%
num_avg = 50:10:100;
% num_avg = [50,100];
EbNodB = 5;
EbNo = 10.^(EbNodB./10);
%%
D_Prob = zeros(length(EbNodB),1);
MD_Prob = zeros(length(EbNodB),1);
FA_Prob = zeros(length(EbNodB),1);
for idx_num_avg = 1:length(num_avg)
    nnum_avg = num_avg(idx_num_avg);
    fprintf('num_avg = %d\n\n',nnum_avg);
    
    p = 0.5*erfc(sqrt(2*EbNo)/2);
    theo_pa = 0.0; theo_rc = 0.0;
    mode_a = 2; %  mode = 1; Unanimous rule, mode = 2; Majority rule
    [Prob_comp min_Prob_usual]=regular_ldpc_analysis_fun3(p,theo_pa,theo_rc,mode_a);
    theo_pa = 1.0; theo_rc = 1.0;
    [Prob_comp max_Prob_usual]=regular_ldpc_analysis_fun3(p,theo_pa,theo_rc,mode_a);
    
    % usual relays
    p2 = max_Prob_usual; % contradiction probability of usual relays
    m2 = nnum_avg*p2; s2 = sqrt(nnum_avg*p2*(1-p2));
    y2 = pdf('normal',0:nnum_avg,m2,s2);
    
    mode = 1;
    rc = 0.15; pa = 1.0;
    n_comp = ceil(nnum_avg*rc);
    n_usual = nnum_avg - n_comp;
    [AveLRContradict] = wb_get_AveLRContradict(rc,pa,EbNodB,nnum_avg,mode);
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
    Knon_Attack_Posi = 1:1:(100*rc);
    Knon_Attack_Posi = Knon_Attack_Posi';
    AttackRelayIndex = est_idx2;
    same = length(intersect(AttackRelayIndex,Knon_Attack_Posi));
    MD = length(setdiff(Knon_Attack_Posi,AttackRelayIndex));
    FA = length(setdiff(AttackRelayIndex,Knon_Attack_Posi));
    
    D_Prob(idx_num_avg) = same/n_comp;
    MD_Prob(idx_num_avg) = MD/n_comp;
    FA_Prob(idx_num_avg) = FA/n_usual;    
end

figure();
plot(num_avg,D_Prob, 'b'); hold on;
plot(num_avg,MD_Prob, 'm');
plot(num_avg,FA_Prob, 'r'); hold off;
legend('Detection Probability','Miss Detection Probability','False Alarm Probability');
% smaller_index = 1:100;
% itr = 10;
% 
% Optimal_T = num_avg +1;
% for i=1:itr
%         
%     %     Termination Test
%     sample_mean = mean(relay_contradict(smaller_index));       
%     sample_var = var(relay_contradict(smaller_index));      
%     binomial_mean = mean(relay_contradict(smaller_index))/num_avg;
%     binomial_var = num_avg*binomial_mean*(1-binomial_mean);         
%     
% %     fprintf('sample mean = %f\n',sample_mean);
% %     fprintf('binomial mean = %f\n\n',binomial_mean);
% %     fprintf('sample variance = %f\n',sample_var);
%     fprintf('binomial variance = %f\n\n',binomial_var);    
%     
%     temp_relay_contradict = relay_contradict(smaller_index);
%     index1 = find(temp_relay_contradict <= sample_mean);
%     index2 = find(temp_relay_contradict > sample_mean);
%     
%     v1 = sum((temp_relay_contradict(index1) - sample_mean).^2);
%     v2 = sum((temp_relay_contradict(index2) - sample_mean).^2);
%             
%     addition (i) = (abs(v1) + abs(v2))/(length(index1)+length(index2));
%     difference (i) = abs(v1 - v2)/(length(index1)+length(index2));
% %     addition (i) = (abs(v1) + abs(v2))/100;
% %     difference (i) = abs(v1 - v2)/100;
%     
%     fprintf('variance 1 = %f\n',v1);
%     fprintf('variacne 2 = %f\n',v2);
%     fprintf('sample variance = %f\n',sample_var);
%     fprintf('difference = %f\n',difference(i));
%     fprintf('addition = %f\n\n',addition(i));
% 
%     n = zeros(1,num_avg+1);
%     for ii=0:num_avg
%         n(ii+1) = histc(relay_contradict(smaller_index),ii);
%     end    
%     figure(); plot(0:1:num_avg,n,'b');hold on;    
%     plot(min_Prob_usual*num_avg,0:0.001:6,'r-');
%     plot(max_Prob_usual*num_avg,0:0.001:6,'c-');
% %     plot(0:1:num_avg,y2,'m-x');
%     
%     n = n / 100;
% 
%     for T = 0:1:num_avg
%         m1 = 0; m2 = 0;
%         for index= 0:1:T
%             m1 = n(index+1)*index + m1 ;
%         end
%         m1 = m1 / sum(n(1:T+1));
%         for index= T+1:1:num_avg
%             m2 = n(index+1)*index + m2 ;
%         end
%         m2 = m2 / sum(n(T+2:num_avg+1));
% 
%         p1= sum(n(1:T+1));
%         p2= 1 - p1;
% 
%         variance(T+1) = p1*p2*(m1-m2)^2;
%     %     variance(T+1) = 0.5*0.5*(m1-m2)^2;
% 
%     end
% 
%     [maximum Optimal_T] = max(abs(variance));
% %     fprintf('maximum between class variacne = %f\n',maximum);
%     
%     
% %     max_var = maximum
%     % Optimal_T = (min(find(abs(variance) == max_var))+max(find(abs(variance) == max_var)) - 1)*0.5/num_avg
%     Optimal_T = (Optimal_T - 1);
% 
%     plot(Optimal_T,0:0.001:6,'g-'); hold off;
% %     xlabel('APC');ylabel('Number of Relays');
% %     legend('APC Histogram of usual realys','APC Histogram of compromised realys','Threshold');
% 
% %     figure(); plot(variance);
% 
%     index_large_real = find(relay_contradict >= Optimal_T)
%     index_small_real = find(relay_contradict < Optimal_T);
%     
%     smaller_index = setdiff(smaller_index, index_large_real);    
%     fprintf('maximum between class variacne = %f\n',maximum);   
%     
% end
% 
% % usual_n = zeros(1,num_avg+1);
% % for i=0:num_avg
% %     usual_n(i+1) = histc(AveLRContradict(100+100*pc+1:200),i);
% % end
% % figure();
% % plot(0:1/num_avg:1,usual_n,'b-o');hold on;
% % axis([0 1 0 6]);
% % 
% % comp_n = zeros(1,num_avg+1);
% % for i=0:num_avg
% %     comp_n(i+1) = histc(AveLRContradict(100+1:100+10),i);
% % end
% % plot(0:1/num_avg:1,comp_n,'r-v');
% % 
% % n = (usual_n + comp_n) / 100;
% % 
% % for T = 0:1:num_avg
% %         
% %     m1 = 0; m2 = 0;
% %     for index= 0:1:T
% %         m1 = n(index+1)*index + m1 ;
% %     end
% %     m1 = m1 / sum(n(1:T+1));
% %     for index= T+1:1:num_avg
% %         m2 = n(index+1)*index + m2 ;
% %     end
% %     m2 = m2 / sum(n(T+2:num_avg+1));
% % 
% %     p1= sum(n(1:T+1));
% %     p2= 1 - p1;
% %     
% %     variance(T+1) = p1*p2*(m1-m2)^2;
% % %     variance(T+1) = 0.5*0.5*(m1-m2)^2;
% %         
% % end
% % 
% % [maximum Optimal_T] = max(abs(variance));
% % 
% % max_var = maximum
% % % Optimal_T = (min(find(abs(variance) == max_var))+max(find(abs(variance) == max_var)) - 1)*0.5/num_avg
% % Optimal_T = (Optimal_T - 1)/num_avg
% % 
% % plot(Optimal_T,0:0.001:6,'r-'); hold off;
% % xlabel('APC');ylabel('Number of Relays');
% % legend('APC Histogram','Threshold');
% % 
% % figure(); plot(variance);
% % 
% % index_large_real = find(relay_contradict/num_avg > Optimal_T)
% % index_small_real = find(relay_contradict/num_avg < Optimal_T);
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis\');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');