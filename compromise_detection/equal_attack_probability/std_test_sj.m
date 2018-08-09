clc;clear;

num_avg = 200;
EbNodB = 5;
mode = 1;
[AveLRContradict] = wb_Main_DAS_NetwrokCoding(0.15,1.0,EbNodB,num_avg,mode);
% load AveLRContradict_5dB AveLRContradict;

n = zeros(1,num_avg+1);
for i=0:num_avg
    n(i+1) = histc(AveLRContradict(101:200),i);
end
n = n/100;
figure();
plot(0:num_avg,n,'x');

relay_contradict = AveLRContradict(101:200);

for T = 0:1:num_avg
        
    m1 = 0; m2 = 0;
    for index= 0:1:T
        m1 = n(index+1)*index + m1 ;
    end
%     m1 = m1 / sum(n(1:T+1));
    for index= T+1:1:num_avg
        m2 = n(index+1)*index + m2 ;
    end
%     m2 = m2 / sum(n(T+2:1000+1));
        
    var1 = 0; var2 = 0;
    for index = 0 : 1 : T
        var1 = n(index+1) * (index - m1)^2 + var1;
    end
%     var1 = var1 / sum(n(1:T+1));
    for index = T+1 :1 : num_avg
        var2 = n(index+1) * (index - m2)^2 + var2;
    end
%     var2 = var2 / sum(n(T+2:1000+1));
    
%     var1 = var1/sum(n(1:T+1)) - m1^2;
%     var2 = var2/sum(n(T+2:1000+1)) - m2^2;     
%     p1= length(find(relay_contradict >= T))/100;
%     p2= length( find(relay_contradict < T))/100;
    
    p1= sum(n(1:T+1));
    p2= 1 - p1;
    
    A = var1 - var2;
    B = 2*(m1*var2 - m2*var1);
    C = var1*m2^2 - var2*m1^2 + 2*var1*var2*log(sqrt(var2)*p1/(sqrt(var1)*p2));
    
    value(T+1) = A*T^2 + B*T + C;
        
end

[minimum Optimal_T] = min(abs(value));

Optimal_T = Optimal_T - 1

index_large_real = find(relay_contradict > Optimal_T)
index_small_real = find(relay_contradict <= Optimal_T)