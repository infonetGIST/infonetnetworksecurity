clc; clear;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis');
%%
% mode = 1; Unanimous rule
% mode = 2; Majority rule
mode = 1;

rc_array = .1:.1:.5;
rc = 0.5;
pa_array = .1:.1:1;
dv = 5; dc = 11;
% dv = 5; dc = 6;

EbNodB = 5; 
EbNo = 10.^(EbNodB./10);

diff_array = zeros(length(EbNodB),length(pa_array));
avg_num_array = zeros(length(EbNodB),length(pa_array));
D_Prob = zeros(length(EbNodB),length(pa_array));
MD_Prob = zeros(length(EbNodB),length(pa_array));
FA_Prob = zeros(length(EbNodB),length(pa_array));

max_num_array = zeros(1, length(EbNodB));

p_ch = zeros(length(EbNodB),1);
for i=1:length(pa_array)
    
    fprintf('Eb/No = %d dB\n',EbNodB);   
    pa = pa_array(i);
    p_ch = 0.5*erfc(sqrt(2*EbNo)/2);
    [pc, pu] = regular_ldpc_analysis_fun3(p_ch,pa,rc,dv,dc,mode);
    
    pc_array(i,1) = pc;
    pu_array(i,1) = pu;    

end

figure();
plot(pa_array,pc_array,'b-x');hold on;
plot(pa_array,pu_array,'r-x');
% plot(pa_array,avg_num_array(3,:),'g-x');
% plot(pa_array,avg_num_array(4,:),'k-x');
% plot(pa_array,avg_num_array(5,:),'m-x');
hold off;
% xlabel('p_{a}');ylabel('Observation Window');
grid on;
% legend('R_c = 0.1','R_c = 0.3','R_c = 0.5','R_c = 0.6','R_c = 0.9');
% legend('R_c = 0.1','R_c = 0.2','R_c = 0.3','R_c = 0.4','R_c = 0.5');
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis');