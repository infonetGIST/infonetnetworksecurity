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

m_large = mean(relay_contradict(index_large));
m_small = mean(relay_contradict(index_small));
m_new = 0.5*(m_large + m_small);

m_temp = m_new;
i = 1;
while(i)    
    
    index_large_temp = find(relay_contradict > m_temp);
    index_small_temp = find(relay_contradict < m_temp);
    
    m_large_temp = mean(relay_contradict(index_large_temp));
    m_small_temp = mean(relay_contradict(index_small_temp));
    m_new_temp = 0.5*(m_large_temp + m_small_temp);    
    
    if abs(m_new_temp - m_temp) < 10^-6
        break;
    end
    m_temp = m_new_temp;
    i = i + 1;
end

m_new_temp
index_large_temp = find(relay_contradict > m_new_temp)
index_small_temp = find(relay_contradict < m_new_temp)

