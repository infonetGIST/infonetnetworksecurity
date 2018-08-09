clc;clear;

num_avg = 1000;
EbNodB = 5;
mode = 1;
pc = 0.11;pa = 0.5;

% [AveLRContradict] = wb_Main_DAS_NetwrokCoding(pc,pa,EbNodB,num_avg,mode);
load AveLRContradict_5dB AveLRContradict;

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

y = (y1*(pc) + y2*(1-pc));

% figure();
plot(0:num_avg,y,'x');
y= n;
for T = 0:1:num_avg    
    
    m1(T+1) = 0; m2(T+1) = 0;
    for index= 0:1:T
        m1(T+1) = y(index+1)*index + m1(T+1) ;
    end
    m1(T+1) = m1(T+1) / sum(y(1:T+1));
    for index= T+1:1:num_avg
        m2(T+1) = y(index+1)*index + m2(T+1) ;
    end
    m2(T+1) = m2(T+1) / sum(y(T+2:num_avg+1));
        
    var1(T+1) = 0; var2(T+1) = 0;
    for index = 0 : 1 : T
        var1(T+1) = y(index+1) * (index - m1(T+1))^2 + var1(T+1);
    end
    var1(T+1) = var1(T+1) / sum(y(1:T+1));
    for index = T+1 :1 : num_avg
        var2(T+1) = y(index+1) * (index - m2(T+1))^2 + var2(T+1);
    end
    var2(T+1) = var2(T+1) / sum(y(T+2:num_avg+1));
    
%     var1 = var1/sum(n(1:T+1)) - m1^2;
%     var2 = var2/sum(n(T+2:1000+1)) - m2^2;     
%     p1= length(find(relay_contradict >= T))/100;
%     p2= length( find(relay_contradict < T))/100;
    
    p1(T+1)= sum(y(1:T+1));
    p2(T+1)= 1 - p1(T+1);
    
    A(T+1) = var1(T+1) - var2(T+1);
    B(T+1) = 2*(m1(T+1)*var2(T+1) - m2(T+1)*var1(T+1));
    C(T+1) = var1(T+1)*m2(T+1)^2 - var2(T+1)*m1(T+1)^2 + 2*var1(T+1)*2*var2(T+1)*log(sqrt(var2(T+1))*p1(T+1)/(sqrt(var1(T+1))*p2(T+1)));    
        
    value(T+1) = abs(A(T+1)*T^2 + B(T+1)*T + C(T+1));
        
end

min_T = ceil(num_avg*Prob_usual);
max_T = ceil(num_avg*Prob_comp);

[minimum Optimal_T] = min(value(min_T:max_T));
Optimal_T = Optimal_T + min_T - 2
plot(Optimal_T, 0:0.0001:0.04);
index_large_real = find(relay_contradict > Optimal_T)
index_small_real = find(relay_contradict <= Optimal_T);

% y = n;
% for T = 0:1:num_avg    
%     
%     m1(T+1) = 0; m2(T+1) = 0;
%     for index= 0:1:T
%         m1(T+1) = y(index+1)*index + m1(T+1) ;
%     end
%     m1(T+1) = m1(T+1) / sum(y(1:T+1));
%     for index= T+1:1:num_avg
%         m2(T+1) = y(index+1)*index + m2(T+1) ;
%     end
%     m2(T+1) = m2(T+1) / sum(y(T+2:1000+1));
%         
%     var1(T+1) = 0; var2(T+1) = 0;
%     for index = 0 : 1 : T
%         var1(T+1) = y(index+1) * (index - m1(T+1))^2 + var1(T+1);
%     end
%     var1(T+1) = var1(T+1) / sum(y(1:T+1));
%     for index = T+1 :1 : num_avg
%         var2(T+1) = y(index+1) * (index - m2(T+1))^2 + var2(T+1);
%     end
%     var2(T+1) = var2(T+1) / sum(y(T+2:1000+1));
%     
% %     var1 = var1/sum(n(1:T+1)) - m1^2;
% %     var2 = var2/sum(n(T+2:1000+1)) - m2^2;     
% %     p1= length(find(relay_contradict >= T))/100;
% %     p2= length( find(relay_contradict < T))/100;
%     
%     p1(T+1)= sum(y(1:T+1));
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
% min_T = ceil(num_avg*Prob_usual);
% max_T = ceil(num_avg*Prob_comp);
% 
% [minimum Optimal_T] = min(value(min_T:max_T));
% Optimal_T = Optimal_T + min_T - 2
% plot(Optimal_T, 0:0.0001:0.04);
% index_large_real = find(relay_contradict > Optimal_T)
% index_small_real = find(relay_contradict <= Optimal_T);


for T = 0:1:num_avg
        
    m1 = 0; m2 = 0;
    for index= 0:1:T
        m1 = n(index+1)*index + m1 ;
    end
    m1 = m1 / sum(n(1:T+1));
    for index= T+1:1:num_avg
        m2 = n(index+1)*index + m2 ;
    end
    m2 = m2 / sum(n(T+2:num_avg+1));

    p1= sum(n(1:T+1));
    p2= 1 - p1;
    
    variance(T+1) = p1*p2*(m1-m2)^2;
%     variance(T+1) = 0.5*0.5*(m1-m2)^2;
        
end

[maximum Optimal_T] = max(abs(variance));
max_var = maximum;
Optimal_range_min = min(find(variance==maximum))-1
Optimal_range_max = max(find(variance==maximum))-1
Optimal_T = Optimal_T - 1;
plot(Optimal_range_min:Optimal_range_max,0.01,'g-');
index_large_real = find(relay_contradict > Optimal_T)
index_small_real = find(relay_contradict <= Optimal_T);

