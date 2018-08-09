%DAS(Distributed Assistant Antenna System) Network performance simulation
clc;clear all;

EEsNO=0:1:5;
% EEsNO = 5;
BERRecord=[];
Ns=100;

for IEsNO=1:length(EEsNO)
    EsNo_dB = EEsNO(IEsNO)
    
    NumOfError = 0;    
    for repeat = 1:1e50        
        % Generate '0' and '1'
        Codeword = randi(1,Ns);
        
        % BPSK
        BPSK_Codeword = 2*Codeword - 1;
        
        % Channel effect
        EsNo = 10^(EsNo_dB/10);
        No = 1/EsNo;
        noise = sqrt(No/2)*(randn(1,Ns));
        yMatrix = BPSK_Codeword + noise;
          
%         received
        y = yMatrix;
            
        y(find(y>0)) = 1;
        y(find(y<0)) = 0;           
        
        NumOfError = length(find(y~=Codeword)) + NumOfError;
        
        %BER breaking criterion
        if (NumOfError >= 1000)
            break;
        end
    end
    
    BER = NumOfError/(repeat*Ns);
    BERRecord = [BERRecord BER];    
    save BERRecord BERRecord
    
    if (BER <= 1e-5)
        break;
    end
    
end %end of EEsN02

figure();
semilogy(EEsNO,BERRecord);
grid on;
% axis([min(EEsNO) max(EEsNO) 1e-6 1e0]);
axis([min(EEsNO) 10 1e-6 1e0]);