% clc;clear;
function [d sum_value] = test(pc,pa,EbNodB,AveLRContradict,num_avg)
% load AveLRContradict_5dB AveLRContradict;


% load Z_comp Z_comp;
% load Z_usual Z_usual;

n = zeros(1,num_avg+1);
for i=0:num_avg
    n(i+1) = histc(AveLRContradict(101:200),i);
end
n = n/100;

% n = n;
% figure();
% plot(0:1000,n,'o');

% Set compromise probabiliy and attack probability
% pc = 0.25;
% pa = 0.2;
% pc = 0.05;
% pa = 1;

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

e = sum(abs(n - y).^2);

% % n = n*100;
% d=0;
% 
% sum_value1 = sum(n(round_m1-round_s1*5:1:round_m1+round_s1*5));
% n(find(n(round_m1-round_s1*5:1:round_m1+round_s1*5)>0)+round_m1-round_s1*5-1) = 0;
% sum_value2 = sum(n(round_m2-round_s2*5:1:round_m2+round_s2*5));
% sum_value = sum_value1 + sum_value2;
% % value1에 포함된거 제외

% if max(AveLRContradict) - round_m1 > 100
%     d = 10^10;
% else
%     for i=round_m1-round_s1*2:1:round_m1+round_s1*2
%         d = n(i)*abs(round_m1-i) + d;
%     end
%     for i=round_m2-round_s2*2:1:round_m2+round_s2*2
%         d = n(i)*abs(round_m2-i) + d;
%     end
% end


% n = n*100;
% d=0;
% sum_value1 = sum(n(round_m1-round_s1*4:1:round_m1+round_s1*4));
% % n(find(n(round_m1-round_s1*4:1:round_m1+round_s1*4)>0)+round_m1-round_s1*4-1) = 0;
% sum_value2 = sum(n(round_m2-round_s2*4:1:round_m2+round_s2*4));
% sum_value = sum_value1 + sum_value2;
% value1에 포함된거 제외

if abs(max(AveLRContradict) - round_m1) > round_s1*6
    d = 10^10;
    sum_value = 10^10;
else
    d = sum(abs(n - y).^2);
    sum_value = sum(n(round_m1-round_s1*4:1:round_m1+round_s1*4))*100;
%     for i=round_m1-round_s1*2:1:round_m1+round_s1*2
%         d = n(i)*abs(round_m1-i) + d;
%     end
%     for i=round_m2-round_s2*2:1:round_m2+round_s2*2
%         d = n(i)*abs(round_m2-i) + d;
%     end
end
