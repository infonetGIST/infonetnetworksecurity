clc;clear;

num_avg = 100 * 10;
EbNodB = 5;
mode = 1;
pc = 0.5;pa = 0.3;

[AveLRContradict] = wb_Main_DAS_NetwrokCoding(pc,pa,EbNodB,num_avg,mode);
% load AveLRContradict_5dB AveLRContradict;

n = zeros(1,num_avg+1);
for i=0:num_avg
    n(i+1) = histc(AveLRContradict(101:200),i);
end
n = n/100;
figure();
plot(0:num_avg,n,'r-x'); hold on;

relay_contradict = AveLRContradict(101:200);

% pc = 0.01;pa = 0.1;
EbNo=10.^(EbNodB./10);
p=0.5*erfc(sqrt(2*EbNo)/2);

[Prob_comp Prob_usual]=regular_ldpc_analysis_fun2(p,pa,pc);

% compromised relays
p1 = Prob_comp; % contradiction probability of compromised relays
m1 = num_avg*p1; s1 = sqrt(num_avg*p1*(1-p1));
y1 = pdf('normal',0:num_avg,m1,s1);
% usual relays
p2 = Prob_usual; % contradiction probability of usual relays
m2 = num_avg*p2; s2 = sqrt(num_avg*p2*(1-p2));
y2 = pdf('normal',0:num_avg,m2,s2);
% Total PDF
y = (y1*(pc) + y2*(1-pc));
% figure();
plot(0:num_avg,y,'x'); hold off;


pc_array = 0:0.01:1;pa_array = 0:0.01:1;

d = ones(length(pc_array),length(pa_array))*10^10;
sum_value = -ones(length(pc_array),length(pa_array));

f = ksdensity(AveLRContradict(101:200),0:num_avg,'function','pdf');
figure();plot(0:num_avg,f);

for i=1:length(pc_array)    
    for j=1:length(pa_array)
        [d(i,j) sum_value(i,j)] = test(pc_array(i),pa_array(j),EbNodB,f,num_avg);
    end
end

[n_row n_col] = size(d);
min_value = min(min(d));
for i=1:n_row
    [A B] = min(d(i,:));
    if A==min_value
        break;
    end
end
est_pc = (i-1)*0.01
% est_pc = sum_value(i,B)*0.01
est_pa = (B-1)*0.01

% [est_Prob_comp est_Prob_usual]=regular_ldpc_analysis_fun2(p,est_pa,est_pc);
% 
% % compromised relays
% est_p1 = est_Prob_comp; % contradiction probability of compromised relays
% est_m1 = num_avg*est_p1; est_s1 = sqrt(num_avg*est_p1*(1-est_p1));
% est_m1_round = round(num_avg*est_p1); est_s1_round = round(sqrt(num_avg*est_p1*(1-est_p1)));
% 
% % est_real_pc = sum(n(est_m1_round-est_s1_round:est_m1_round+est_s1_round))
% % 
% % [est_Prob_comp est_Prob_usual]=regular_ldpc_analysis_fun2(p,est_pa,est_real_pc);
% % % compromised relays
% % est_p1 = est_Prob_comp; % contradiction probability of compromised relays
% % est_m1 = num_avg*est_p1; est_s1 = sqrt(num_avg*est_p1*(1-est_p1));
% % est_m1_round = round(num_avg*est_p1); est_s1_round = round(sqrt(num_avg*est_p1*(1-est_p1)));
% 
% est_y1 = pdf('normal',0:num_avg,est_m1,est_s1);
% % usual relays
% est_p2 = est_Prob_usual; % contradiction probability of usual relays
% est_m2 = num_avg*est_p2; est_s2 = sqrt(num_avg*est_p2*(1-est_p2));
% est_y2 = pdf('normal',0:num_avg,est_m2,est_s2);
% 
% est_y = (est_y1*(est_pc) + est_y2*(1-est_pc));
% 
% figure();
% plot(0:num_avg,n,'r-x'); hold on;
% plot(0:num_avg,est_y,'x');
% 
% for T = 0:1:num_avg    
%     
%     m1(T+1) = 0; m2(T+1) = 0;
%     for index= 0:1:T
%         m1(T+1) = est_y(index+1)*index + m1(T+1) ;
%     end
%     m1(T+1) = m1(T+1) / sum(est_y(1:T+1));
%     for index= T+1:1:num_avg
%         m2(T+1) = est_y(index+1)*index + m2(T+1) ;
%     end
%     m2(T+1) = m2(T+1) / sum(est_y(T+2:num_avg+1));
%         
%     var1(T+1) = 0; var2(T+1) = 0;
%     for index = 0 : 1 : T
%         var1(T+1) = est_y(index+1) * (index - m1(T+1))^2 + var1(T+1);
%     end
%     var1(T+1) = var1(T+1) / sum(est_y(1:T+1));
%     for index = T+1 :1 : num_avg
%         var2(T+1) = est_y(index+1) * (index - m2(T+1))^2 + var2(T+1);
%     end
%     var2(T+1) = var2(T+1) / sum(est_y(T+2:num_avg+1));
%         
%     p1(T+1)= sum(est_y(1:T+1));
%     p2(T+1)= 1 - p1(T+1);
%     
%     A(T+1) = var1(T+1) - var2(T+1);
%     B(T+1) = 2*(m1(T+1)*var2(T+1) - m2(T+1)*var1(T+1));
%     C(T+1) = var1(T+1)*m2(T+1)^2 - var2(T+1)*m1(T+1)^2 + 2*var1(T+1)*2*var2(T+1)*log(sqrt(var2(T+1))*p1(T+1)/(sqrt(var1(T+1))*p2(T+1)));    
%         
%     value(T+1) = abs(A(T+1)*T^2 + B(T+1)*T + C(T+1));
%         
% end
% 
% % min_T = ceil(num_avg*est_Prob_usual);
% % max_T = floor(num_avg*est_Prob_comp);
% % min_T = min(find(n~=0));
% % max_T = max(find(n~=0));
% 
% min_T = max(ceil(num_avg*est_Prob_usual), min(find(n~=0)));
% max_T = min(floor(num_avg*est_Prob_comp), max(find(n~=0)));
% 
% [minimum Optimal_T] = min(value(min_T:max_T));
% Optimal_T = Optimal_T + min_T - 2
% % plot(Optimal_T, 0:0.0001:0.03,'m-x');
% index_large_real = find(relay_contradict > Optimal_T)
% index_small_real = find(relay_contradict <= Optimal_T);
% 
% % y = n;
% % for T = 0:1:num_avg    
% %     
% %     m1(T+1) = 0; m2(T+1) = 0;
% %     for index= 0:1:T
% %         m1(T+1) = y(index+1)*index + m1(T+1) ;
% %     end
% %     m1(T+1) = m1(T+1) / sum(y(1:T+1));
% %     for index= T+1:1:num_avg
% %         m2(T+1) = y(index+1)*index + m2(T+1) ;
% %     end
% %     m2(T+1) = m2(T+1) / sum(y(T+2:num_avg+1));
% %         
% %     var1(T+1) = 0; var2(T+1) = 0;
% %     for index = 0 : 1 : T
% %         var1(T+1) = y(index+1) * (index - m1(T+1))^2 + var1(T+1);
% %     end
% %     var1(T+1) = var1(T+1) / sum(y(1:T+1));
% %     for index = T+1 :1 : num_avg
% %         var2(T+1) = y(index+1) * (index - m2(T+1))^2 + var2(T+1);
% %     end
% %     var2(T+1) = var2(T+1) / sum(y(T+2:num_avg+1));
% %     
% % %     var1 = var1/sum(n(1:T+1)) - m1^2;
% % %     var2 = var2/sum(n(T+2:1000+1)) - m2^2;     
% % %     p1= length(find(relay_contradict >= T))/100;
% % %     p2= length( find(relay_contradict < T))/100;
% %     
% %     p1(T+1)= sum(y(1:T+1));
% %     p2(T+1)= 1 - p1(T+1);
% %     
% %     A(T+1) = var1(T+1) - var2(T+1);
% %     B(T+1) = 2*(m1(T+1)*var2(T+1) - m2(T+1)*var1(T+1));
% %     C(T+1) = var1(T+1)*m2(T+1)^2 - var2(T+1)*m1(T+1)^2 + 2*var1(T+1)*2*var2(T+1)*log(sqrt(var2(T+1))*p1(T+1)/(sqrt(var1(T+1))*p2(T+1)));    
% %         
% %     value(T+1) = abs(A(T+1)*T^2 + B(T+1)*T + C(T+1));
% %         
% % end
% % 
% % % min_T = min(find(n~=0));
% % % max_T = max(find(n~=0));
% % min_T = max(ceil(num_avg*est_Prob_usual), min(find(n~=0)));
% % max_T = min(floor(num_avg*est_Prob_comp), max(find(n~=0)));
% % 
% % [minimum Optimal_T] = min(value(min_T:max_T));
% % Optimal_T = Optimal_T + min_T - 2
% % plot(Optimal_T, 0:0.0001:0.04,'g-x');
% % index_large_real = find(relay_contradict > Optimal_T)
% % index_small_real = find(relay_contradict <= Optimal_T);
% 
% 
% 
% % for T = 0:1:num_avg
% %         
% %     m1 = 0; m2 = 0;
% %     for index= 0:1:T
% %         m1 = est_y(index+1)*index + m1 ;
% %     end
% %     m1 = m1 / sum(n(1:T+1));
% %     for index= T+1:1:num_avg
% %         m2 = est_y(index+1)*index + m2 ;
% %     end
% %     m2 = m2 / sum(est_y(T+2:num_avg+1));
% % 
% %     p1= sum(est_y(1:T+1));
% %     p2= 1 - p1;
% %     
% %     variance(T+1) = p1*p2*(m1-m2)^2;       
% % end
% % max_var = max(variance(find(variance ~= inf)));
% % Optimal_T = find(variance == max_var) - 1
% % plot(Optimal_T,0:.00001:0.06,'g-');
% % index_large_real = find(relay_contradict > Optimal_T)
% % index_small_real = find(relay_contradict <= Optimal_T);
