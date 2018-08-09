clc;clear;
%%

% EEsNO = -4:1:9;
EEsNO = 10;
% EEsNO = 7;
BERRecord = [];
Ns = 100;

for IEsNO=1:length(EEsNO)
    EsNodB = EEsNO(IEsNO);
    fprintf('Es/No = %d [dB]\n',EsNodB);
    
    NumOfError = 0;    
    for repeat = 1:1e12        
        % Generate '0' and '1'
        Codeword = randsrc(1, Ns, [0 1]);
        
        % BPSK
        BPSK_Codeword = 2*Codeword - 1;
        
        % Channel effect
        EsNo = 10^(EsNodB/10);
        No = 1/EsNo;
        noise = sqrt(No/2)*(randn(1,Ns));
        yMatrix = BPSK_Codeword + noise;
        
%         yMatrix = awgn(BPSK_Codeword, EsNodB);
        
%         received
        y = yMatrix;
            
        y(find(y>0)) = 1;
        y(find(y<0)) = 0;
            
        NumOfError = sum(mod(y + Codeword,2)) + NumOfError;
%         NumOfError = length(find(y~=Codeword)) + NumOfError;
        
        %BER breaking criterion
        if (NumOfError >= 1000)
            break;
        end
    end
    
    BER = NumOfError/(repeat*Ns);
    BERRecord = [BERRecord BER];    
%     save BERRecord BERRecord
    BER_th(IEsNO,:) = 1/2*erfc(sqrt(EsNo));
    if (BER <= 1e-4)
        break;
    end
    
end %end of EEsN02



figure();
semilogy(EEsNO,BERRecord); hold on;
semilogy(EEsNO,BER_th); hold off;
grid on;
% axis([min(EEsNO) max(EEsNO) 1e-6 1e0]);
axis([0 10 1e-6 1e0]);
