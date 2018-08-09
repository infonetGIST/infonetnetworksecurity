% clc;clear;
function [d sum_value] = test(pc,pa,EbNodB,n,num_avg)

% n = zeros(1,num_avg+1);
% for i=0:num_avg
%     n(i+1) = histc(AveLRContradict(101:200),i);
% end
% n = n/100;



% Obtain theoretical contradiction probability

% EbNodB = 5;
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

% multiplying prior probability

y = (y1*(pc) + y2*(1-pc));
% r_y = ceil(y*100);
% y = r_y/100;
% figure();
% plot(0:1000,y,'x');
round_m1 = round(m1)+1; round_s1 = round(s1);
round_m2 = round(m2)+1; round_s2 = round(s2);


% % using mean difference
% max_nonzero = max(find(n~=0));
% i = max_nonzero;
% while i
%     sum_value = sum(n(i:num_avg));
%     if sum_value >= pc
%         break;
%     end
%     i = i -1;
% end
% temp_index1 = i:max_nonzero;
% v1 = sum(abs(temp_index1-1 - round_m1).^2.*n(temp_index1));
% temp_index2 = 1:i-1;
% v2 = sum(abs(temp_index2-1 - round_m2).^2.*n(temp_index2));
% d = v1 + v2;
% 
% if abs(max(find(n~=0)) - round_m1) > 5*round_s1
%     d = 10^10;
% end


% using square error

temp_n = n;
sum_value1 = sum(temp_n(round_m1-round_s1*3:1:round_m1+round_s1*3));
temp_n(find(temp_n(round_m1-round_s1*3:1:round_m1+round_s1*3)>0) + round_m1-round_s1*3-1) = 0;
sum_value2 = sum(temp_n(round_m2-round_s2*3:1:round_m2+round_s2*3));

if sum_value1 <= pc*0.65 || sum_value2 <= (1-pc)*0.65
% if sum_value1 <= pc*0.65 || sum_value2 <= (1-pc)*0.65 || abs(max(find(n~=0)) - round_m1) > 6*round_s1    
    d = 10^10;
    sum_value = 10^10;
else
    d = sum(abs(n - y).^2);
    sum_value = sum(n(round_m1-round_s1*4:1:round_m1+round_s1*4))*100;

end
