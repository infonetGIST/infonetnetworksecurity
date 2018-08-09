clc; clear;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/contradiction_analysis');
%%
% EbNodB = 5;
% EbNo = 10.^(EbNodB./10);
% p = 0.5*erfc(sqrt(2*EbNo)/2);
% pa = 0:0.01:1;
% compromise rate
% rc = 0.15;

% Ns = 50; Nr = 100; dv = 10; dc = 6;
% Ns = 100; Nr = 100; dv = 8; dc = 9;
Ns = 200; Nr = 100; dv = 6; dc = 13;
%%
code_rate = Ns / (Ns+Nr);
% EbNodB = 8;
% EbNo = 10.^(EbNodB./10);
% EsNo = EbNo * code_rate;
% EsNodB = 10*log10(EsNo);
EsNodB = 5;
EsNo = 10.^(EsNodB./10);
%%
% mode = 1; Unanimous rule, mode = 2; Majority rule
mode = 1;

rc_array = .1:.1:.5;
pa_array = .1:.1:1;
% rc_array = .5;
% pa_array = 1;

% EbNodB = 5;
% EbNodB = 7;
% EbNo = 10.^(EbNodB./10);

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
%     p_ch(index) = 0.5*erfc(sqrt(2*EbNo(index))/2);
    p_ch(index) = 0.5*erfc(sqrt(EsNo(index)));

%     Rayleigh fading channel
%     tmp = 0;
%     for k =1:1000
%         rho = sqrt(0.5)*(randn(1,1) + 1j*randn(1,1));
%         tmp = 0.5* erfc(abs(rho)*sqrt(EbNo(index))) + tmp;
%     end
%     p_ch(index) = tmp/k;    
    
    for a = 1:length(rc_array)
        rc = rc_array(a);    
        for b = 1:length(pa_array)
            pa = pa_array(b);
            pr = p_ch + pa - 2*p_ch*pa;

            if EsNodB == 3
                num_array = 100:10:1000000;
            elseif EsNodB == 4
                num_array = 1:1:4000;
            elseif EsNodB == 5
                num_array = 10:1:3000;
%                 num_array = 1000;
            elseif EsNodB == 6
                num_array = 10:1:3000;
            else
                num_array = 10:1:100^4;
            end
            
            % finding threshold to make false alarm probability less than alpha            
            itr = 1;
            while itr
                num_avg = num_array(itr);
                [pc, pu] = regular_ldpc_analysis_fun3(p_ch(index),pa,rc,dv,dc,mode);
                
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
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\contradiction_analysis');
% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/contradiction_analysis');
