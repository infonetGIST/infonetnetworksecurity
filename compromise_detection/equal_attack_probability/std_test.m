clc;clear;

num_avg = 1000;
EbNodB = 5;
% [AveLRContradict] = wb_Main_DAS_NetwrokCoding(0.01,0.1,EbNodB,num_avg);
load AveLRContradict_5dB AveLRContradict;

n = zeros(1,num_avg+1);
for i=0:num_avg
    n(i+1) = histc(AveLRContradict(101:200),i);
end
n = n/100;
figure();
plot(0:num_avg,n,'x');

relay_contradict = AveLRContradict(101:200);
m = mean(relay_contradict);
index_large = find(relay_contradict > m);
index_small = find(relay_contradict < m);
std_large = std(relay_contradict(index_large));
std_small = std(relay_contradict(index_small));
std_sum = std_large + std_small;

m_new = 0.5*(mean(relay_contradict(index_large)) +mean(relay_contradict(index_small)));
m_old = m_new;
index_large_min = index_large;
index_small_min = index_small;
std_large_min = std_large;
std_small_min = std_small;
std_sum_min = std_sum;

i = 1;
while(i)
    
    m_temp = m_new;
    index_large_temp = find(relay_contradict > m_temp);
    index_small_temp = find(relay_contradict < m_temp);
    std_large_temp = std(relay_contradict(index_large_temp));
    std_small_temp = std(relay_contradict(index_small_temp));
    std_sum_temp = std_large_temp + std_small_temp;
    
%   move to smaller standard deviation  
    if std_sum_temp < std_sum_min
        index_large_min = index_large_temp;
        index_small_min = index_small_temp;
        std_large_min = std_large_temp;
        std_small_min = std_small_temp;
        std_sum_min = std_sum_temp;
    else
%         check next standard deviation is smaller or not
        m_new = 0.5*(mean(relay_contradict(index_large_temp)) + mean(relay_contradict(index_small_temp)));
        index_large_new = find(relay_contradict > m_new);
        index_small_new = find(relay_contradict < m_new);
        std_large_new = std(relay_contradict(index_large_new));
        std_small_new = std(relay_contradict(index_small_new));
        std_sum_new = std_large_new + std_small_new;        
        
        if std_sum_new < std_sum_min
            index_large_min = index_large_new;
            index_small_min = index_small_new;
            std_large_min = std_large_new;
            std_small_min = std_small_new;
            std_sum_min = std_sum_new;
            i = i + 1;
            continue;
        else
            m_forward = m_old;
            std_sum_forward = std_sum_min;
            index_large_forward = index_large_min;
            index_small_forward = index_small_min;
            break;
        end
    end
    
    m_old = m_temp;
    m_new = 0.5*(mean(relay_contradict(index_large_min)) +mean(relay_contradict(index_small_min)));    
    if m_new > 1000
        m_forward = m_old;
        std_sum_forward = std_sum_min;
        index_large_forward = index_large_min;
        index_small_forward = index_small_min;
        break;
    end
    i = i + 1;
end

