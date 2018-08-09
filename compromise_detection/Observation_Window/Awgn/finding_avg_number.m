clc; clear;

% rc_array = 0.1:0.1:0.5;
rc_array = [0.1 0.3 0.5 0.7 0.9];
pa_array = 0.1:0.1:1;

EsN0_dB = 5;
EsN0 = 10.^(EsN0_dB./10);
EbN0 = 2*EsN0;
EbN0_dB = 10*log10(EbN0);

diff_array = zeros(length(EsN0_dB),length(pa_array));
avg_num_array = zeros(length(EsN0_dB),length(pa_array));
D_Prob = zeros(length(EsN0_dB),length(pa_array));
MD_Prob = zeros(length(EsN0_dB),length(pa_array));
FA_Prob = zeros(length(EsN0_dB),length(pa_array));

max_num_array = zeros(1, length(EsN0_dB));

p_ch = zeros(length(EsN0_dB),1);
for index=1:length(EsN0_dB)
    
    fprintf('Es/No = %d dB\n',EsN0_dB);
    
%     tmp = 0;
%     for k =1:1000
%         rho = sqrt(0.5)*(randn(1,1) + 1j*randn(1,1));
%         tmp = 0.5* erfc(abs(rho)*sqrt(EsN0(index))) + tmp;
%     end
%     p_ch(index) = tmp/k;
%     Only AWGN Channel
    p_ch(index) = 0.5*erfc(sqrt(EsN0(index)));

    for a = 1:length(rc_array)
        rc = rc_array(a);    
        for b = 1:length(pa_array)
            pa = pa_array(b);
            p_r = p_ch + pa - 2*p_ch*pa;

            if EsN0_dB == 3
                num_array = 100:10:1000000;
            elseif EsN0_dB == 4
                num_array = 1:1:4000;
            elseif EsN0_dB ==5
                num_array = 1:1:300000;
            else
                num_array = 100:10:100^4;
            end            
            
            % finding threshold to make false alarm probability less than alpha            
            itr = 1;
            while itr
                num_avg = num_array(itr);
                [pc pu]=regular_ldpc_analysis_fun2(p_ch(index),p_r,rc);
%                 pc = 0.3; pu = 0.2;
                
                x = 0:num_avg;
%                 PDF of normal realys
                y1 = binopdf(x,num_avg,pu)*(1-rc);                

                s = 0; i = max(x) + 1;
                while i    
                    s = y1(i) + s;    
                    if s >= 10^-4
                        break;
                    end
                    i = i - 1;
                end
                p_F = sum(y1(i:num_avg+1))/(1-rc);
                
                y2 = binopdf(x,num_avg,pc)*rc;                
                p_D = sum(y2(i:num_avg+1))/rc;                


                if abs(p_D - 1) < 10^-4 || num_array(itr)==10^4
%                 if num_array(itr)==10^4
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
plot(pa_array,avg_num_array(5,:),'m-x'); hold off;
xlabel('p_{a}');ylabel('Observation Window');
grid on;
legend('R_c = 0.1','R_c = 0.3','R_c = 0.5','R_c = 0.6','R_c = 0.9');