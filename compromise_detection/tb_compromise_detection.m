clc;clear;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis\');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
%%
EbNodB = 5;

EbNo = 10.^(EbNodB./10);
p = 0.5*erfc(sqrt(2*EbNo)/2);
theo_pa = 0.0; theo_rc = 0.0;
mode_a = 2; %  mode = 1; Unanimous rule, mode = 2; Majority rule
[Prob_comp min_Prob_usual]=regular_ldpc_analysis_fun3(p,theo_pa,theo_rc,mode_a);
theo_pa = 1.0; theo_rc = 1.0;
[Prob_comp max_Prob_usual]=regular_ldpc_analysis_fun3(p,theo_pa,theo_rc,mode_a);

num_avg = 100;
% usual relays
p2 = max_Prob_usual; % contradiction probability of usual relays
m2 = num_avg*p2; s2 = sqrt(num_avg*p2*(1-p2));
y2 = pdf('normal',0:num_avg,m2,s2);
%%
mode = 1;
rc = 0.12; pa = 1.0;
n_comp = num_avg*rc;
n_usual = num_avg - n_comp;
% [AveLRContradict] = wb_get_AveLRContradict(rc,pa,EbNodB,num_avg,mode);
load AveLRContradict_5dB AveLRContradict;

relay_contradict = AveLRContradict(101:200)';
%% K-means clustering
ini_point = [0;100];
opts = statset('Display','final','MaxIter',100);
[idx,C,sumd,D] = kmeans(relay_contradict,2,'Distance','cityblock',...
'Replicates',1,'Start',ini_point,'Options',opts);
%     'Replicates',1,'Options',opts);    
est_idx1 = find(idx==1);
est_idx2 = find(idx==2);

figure;
plot(est_idx1,relay_contradict(est_idx1,1),'r.','MarkerSize',12); hold on;
plot(est_idx2,relay_contradict(est_idx2,1),'b.','MarkerSize',12); hold off;
%%
Knon_Attack_Posi = 1:1:n_comp;
AttackRelayIndex = est_idx2;
temp=0;
for i=1:1:length(Knon_Attack_Posi)
    for j=1:1:length(AttackRelayIndex)
        if Knon_Attack_Posi(i) == AttackRelayIndex(j)
            temp = temp +1;
        end
    end
    same = temp;
    MD = length(Knon_Attack_Posi) - temp;
    FA = length(AttackRelayIndex) - temp;
end
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis\');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');