clc;clear;

num_avg = 1000;
EbNodB = 5;
mode = 2;
% [AveLRContradict] = wb_Main_DAS_NetwrokCoding(0.5,0.5,EbNodB,num_avg,mode);
load AveLRContradict_5dB AveLRContradict;

EbNo=10.^(EbNodB./10);
p=0.5*erfc(sqrt(2*EbNo)/2);
pa = 1.0; pc = 1.0;
[Prob_comp Prob_usual]=regular_ldpc_analysis_fun2(p,pa,pc);

n = zeros(1,num_avg+1);
for i=0:num_avg
    n(i+1) = histc(AveLRContradict(101:200),i);
end
n = n/100;
figure();
plot(0:num_avg,n,'x');hold on;
plot(Prob_usual*num_avg,0:0.0001:0.05);hold off;

relay_contradict = AveLRContradict(101:200);
m = mean(relay_contradict);
index_large = find(relay_contradict > m);
index_small = find(relay_contradict < m);
std_large = std(relay_contradict(index_large));
std_small = std(relay_contradict(index_small));
std_sum = std_large + std_small;

% m_new = 0.5*(mean(relay_contradict(index_large)) +mean(relay_contradict(index_small)));
m_old = m;
m_min = m_old;
index_large_min = index_large;
index_small_min = index_small;
std_large_min = std_large;
std_small_min = std_small;
std_sum_min = std_sum;

i = 1;
while(i)    
    m_temp = m_old + 1;
    index_large_temp = find(relay_contradict > m_temp);
    index_small_temp = find(relay_contradict < m_temp);
    std_large_temp = std(relay_contradict(index_large_temp));
    std_small_temp = std(relay_contradict(index_small_temp));
    std_sum_temp = std_large_temp + std_small_temp;
    
%   move to smaller standard deviation  
    if std_sum_temp < std_sum_min
        m_min = m_temp;
        index_large_min = index_large_temp;
        index_small_min = index_small_temp;
        std_large_min = std_large_temp;
        std_small_min = std_small_temp;
        std_sum_min = std_sum_temp;
    end
    
    m_old = m_temp;
    
%     if m_old >= 1000 || isempty(index_large_temp) || isempty(index_small_temp)
    if m_old >= 950
        m_forward = m_min;
        std_sum_forward = std_sum_min;
        index_large_forward = index_large_min;
        index_small_forward = index_small_min;
        break;
    end
    i = i + 1;
end

% m_new = 0.5*(mean(relay_contradict(index_large)) +mean(relay_contradict(index_small)));
m_old = m;
index_large_min = index_large;
index_small_min = index_small;
std_large_min = std_large;
std_small_min = std_small;
std_sum_min = std_sum;
j = 1;
while(j)    
    m_temp = m_old - 1;
    index_large_temp = find(relay_contradict > m_temp);
    index_small_temp = find(relay_contradict < m_temp);
    std_large_temp = std(relay_contradict(index_large_temp));
    std_small_temp = std(relay_contradict(index_small_temp));
    std_sum_temp = std_large_temp + std_small_temp;
    
%   move to smaller standard deviation  
    if std_sum_temp < std_sum_min
        m_min = m_temp;
        index_large_min = index_large_temp;
        index_small_min = index_small_temp;
        std_large_min = std_large_temp;
        std_small_min = std_small_temp;
        std_sum_min = std_sum_temp;
    end
    
    m_old = m_temp;
    
%     if m_old >= 1000 || isempty(index_large_temp) || isempty(index_small_temp)
    if m_old <= 10
        m_backward = m_min;
        std_sum_backward = std_sum_min;
        index_large_backward = index_large_min;
        index_small_backward = index_small_min;
        break;
    end
    j = j + 1;
end

if std_sum_forward < std_sum_backward
    est_m = m_forward;
    index_large_est = index_large_forward;
    index_small_est = index_small_forward;
else
    est_m = m_backward;
    index_large_est = index_large_backward;
    index_small_est = index_small_backward;
end
% if abs(std_sum_forward - std_sum_backward) < 1 && m_backward > 81.7
%     index_large_est = 1:100;
%     index_small_est = 0;
% end
% if abs(std_sum_forward - std_sum_backward) < 1 && m_backward < 81.7
%     index_large_est = 0;
%     index_small_est = 1:100;
% end
est_m
index_large_est
index_small_est