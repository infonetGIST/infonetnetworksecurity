clc;clear;
% Rc = 0.05
% num_avg = 840;
% Rc = 0.15
num_avg = 906;
EbNodB = 5;
mode = 1;
pc = 0.15;pa = 1.0;
[AveLRContradict] = wb_Main_DAS_NetwrokCoding(pc,pa,EbNodB,num_avg,mode);
% load AveLRContradict_5dB AveLRContradict;

EbNo=10.^(EbNodB./10);
p=0.5*erfc(sqrt(2*EbNo)/2);
% pa = 1.0; pc = 1.0;
% [Prob_comp Prob_usual]=regular_ldpc_analysis_fun2(p,pa,pc);

n = zeros(1,num_avg+1);
for i=0:num_avg
    n(i+1) = histc(AveLRContradict(101:200),i);
end
% n = n/100;
figure();
% plot(0:1/num_avg:1,n,'x');hold on;
% plot(Prob_usual*num_avg,0:0.0001:0.05);hold off;
plot(0:1/num_avg:1,n,'b');hold on;
% bar(0:1/num_avg:1,n,'b');hold on;
axis([0 1 0 6]);


% relay_contradict = AveLRContradict(101:200);
% 
% n = n / 100;
% for T = 0:1:num_avg
%         
%     m1 = 0; m2 = 0;
%     for index= 0:1:T
%         m1 = n(index+1)*index + m1 ;
%     end
%     m1 = m1 / sum(n(1:T+1));
%     for index= T+1:1:num_avg
%         m2 = n(index+1)*index + m2 ;
%     end
%     m2 = m2 / sum(n(T+2:num_avg+1));
% 
%     p1= sum(n(1:T+1));
%     p2= 1 - p1;
%     
%     variance(T+1) = p1*p2*(m1-m2)^2;
% %     variance(T+1) = 0.5*0.5*(m1-m2)^2;
%         
% end
% 
% [maximum Optimal_T] = max(abs(variance));
% 
% max_var = maximum
% Optimal_T = (min(find(abs(variance) == max_var))+max(find(abs(variance) == max_var)) - 1)*0.5/num_avg
% % Optimal_T = (Optimal_T - 1)/num_avg
% 
% 
% plot(Optimal_T,0:0.001:6,'r-');
% xlabel('APC');ylabel('Number of Relays');
% legend('APC Histogram','Threshold');
% 
% index_large_real = find(relay_contradict/num_avg > Optimal_T)
% % index_small_real = find(relay_contradict <= Optimal_T)
% % % 
% % % [a b]=graythresh(n);
% % % b=b*100
% % % index_large_real = find(relay_contradict > b)
% % index_small_real = find(relay_contradict <= b)
% % 
