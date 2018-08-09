% ***************************************************************** 
% COPYRIGHT (c) 2018 Heung-No Lee, and Woong-Bi Lee. 
% E-mail: heungno@gist.ac.kr, woongbi.lee@gmail.com
% Affiliation: INFONET Laboratory, Gwangju Institute of Science and
% Technology (GIST), Republic of Korea
% homepage: http://infonet.gist.ac.kr
% *****************************************************************  
% filename: Observation_window.m
% this script generates observation window to satisfy the false alarm probability
% *****************************************************************
%% Parameters
% Ns: # of sources, Nr: # of relays, dv: # of edges per variable node, dc: # of edges per check node
% case 1
% Ns = 50; Nr = 100; dv = 10; dc = 6; 
% case 2
Ns = 100; Nr = 100; dv = 8; dc = 9;
% case 3
% Ns = 200; Nr = 100; dv = 6; dc = 13; 
EsNodB = 5;
rc_array = .1:.1:.5; % compromise rate
pa_array = .1:.1:1; % attack probability
mode = 1; % mode = 1; Unanimity rule, mode = 2; Majority rule
%%
EsNo = 10.^(EsNodB./10);

diff_array = zeros(length(EsNodB),length(pa_array));
avg_num_array = zeros(length(EsNodB),length(pa_array));
D_Prob = zeros(length(EsNodB),length(pa_array));
MD_Prob = zeros(length(EsNodB),length(pa_array));
FA_Prob = zeros(length(EsNodB),length(pa_array));

max_num_array = zeros(1, length(EsNodB));

p_ch = zeros(length(EsNodB),1);
for index=1:length(EsNodB)
    
    fprintf('Eb/No = %d dB\n',EsNodB);   

%   AWGN Channel
    p_ch(index) = 0.5*erfc(sqrt(EsNo(index)));
    
    for a = 1:length(rc_array)
        rc = rc_array(a);    
        for b = 1:length(pa_array)
            pa = pa_array(b);
            pr = p_ch + pa - 2*p_ch*pa;

            num_array = 10:1:3000;
            
            % finding threshold to make false alarm probability less than alpha            
            itr = 1;
            while itr
                num_avg = num_array(itr);
                [pc, pu] = regular_ldpc_analysis_fun2(p_ch(index),pa,rc,dv,dc,mode);
                
                x = 0:num_avg;
                y1 = binopdf(x,num_avg,pu);

                i = length(y1);
                while i
                    p_F = sum(y1(i:end));
                    if p_F >= 1e-3
                        break;
                    end
                    i = i - 1;
                end
                
                y2 = binopdf(x,num_avg,pc);                
                p_D = sum(y2(i:end));

                if p_D > 1 - 1e-3 || num_array(itr)==10^4
                    break;
                end
                itr = itr + 1;
            end
            avg_num_array(a, b) = num_avg;
        end
    end    
    max_num_array(index) = max(max(avg_num_array));    
end

figure();
plot(pa_array,avg_num_array(1,:),'b-x');hold on;
plot(pa_array,avg_num_array(2,:),'r-x');
plot(pa_array,avg_num_array(3,:),'g-x');
plot(pa_array,avg_num_array(4,:),'k-x');
plot(pa_array,avg_num_array(5,:),'m-x');
hold off;
xlabel('p_{a}');ylabel('Observation Window');
grid on;
% legend('R_c = 0.1','R_c = 0.3','R_c = 0.5','R_c = 0.6','R_c = 0.9');
legend('R_c = 0.1','R_c = 0.2','R_c = 0.3','R_c = 0.4','R_c = 0.5');