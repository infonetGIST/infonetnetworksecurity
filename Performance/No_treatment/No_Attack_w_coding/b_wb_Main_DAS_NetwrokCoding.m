%DAS(Distributed Assistant Antenna System) Network performance simulation
clc;clear;close all;
%%
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');

% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/LDPC_decoder');
%%
% case 1
% Ns = 50; Nr = 100; dv = 10; dc = 6;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');
% case 2
% Ns = 100; Nr = 100; dv = 8; dc = 9;
Ns = 100; Nr = 100; dv = 4; dc = 5;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns100_icr8_dv8_dc9');
% case 3
% Ns = 200; Nr = 100; dv = 6; dc = 13;
% % addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns200_icr12_dv6_dc13');
code_rate = Ns / (Ns+Nr);
%%
% Ns=100; %# of Senders
% Nd=100; %# of distributed relays
% ICR=5;%# of incoming connection per relay 
% %Enter 1 for no time-domain coding
% n=1;  %# of output coded bits in time-domain
k=1;  %# of input bits in time-domain
% H_j=1;%# of 1's in a column in the H matrix
% H_k=1;%# of 1's in a row i the H matrix
%%%%# of inner iterations in LDPC decoder
Num_iteration=20; 

AttackPercent=0; %0.05 %# % of relay nodes are attached
Pinsist=1;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2

EbNodB_array = -2:1:5;
% EbNodB_array = -2:1:5;

BERRecord=[];
AveLRContradict=[];

for idx=1:length(EbNodB_array)
%     EbNodB = EbNodB_array(idx); %unit in dB
%     EbNo = 10.^(EbNodB./10);
%     EsNo = EbNo * code_rate;
%     EsNodB = 10*log10(EsNo);    
%     No = 1/EsNo;
    
    EbNodB = EbNodB_array(idx); %unit in dB
    EbNo = 10.^(EbNodB./10);
    RNo = 1/EbNo;
        
    fprintf('Eb/No = %d [dB]\n',EbNodB);
    fprintf('Noise power = %f \n',RNo);

    NumOfError = 0;    
    noAveSamples = 100; %# of test samples
    
    for repeat = 1:1e6
        MessageMatrix = randsrc(noAveSamples, Ns, [0 1]);
        LRContradict = 0;
        
        for NumOfTx = 1:noAveSamples
            Codeword = MessageMatrix(NumOfTx,:);
            MeshCodeword=[];
            
            for Inc = 1:size(Codeword, 1) % Inc : 1 ~ NumOfTx
                MeshCodewordTemp = mod(Codeword(Inc,:)*G_NC ,2); % Generating Mesh Codeword
                MeshCodeword = [MeshCodeword; MeshCodewordTemp];
            end
            
            %Channel effect
%             EsN02 = EEsNO2(IEsNO2); %unit in dB
%             yMatrix = awgn((2*MeshCodeword - 1), EsNO2);
            noise = sqrt(RNo/2)*randn(1,Ns+Nr);
            signal = 2*MeshCodeword - 1;
            yMatrix = signal + noise;
            
            %LDPC decodor            
            y=reshape(yMatrix', 1, size(yMatrix,1)*size(yMatrix,2));
            
%             EsNO2_lin = 10^(EsNO2/10);
%             LR_f = (2*EsNO2_lin).*y;
            LR_f = (2*EbNo).*y;
%             LR_f = (4*EbNo).*y;
            LR_fMatrix(NumOfTx,:) = LR_f;   
        end

        for index = 1:size(LR_fMatrix, 1)            
            LR_f = LR_fMatrix(index,:);
            %decoding process
            LR_f((LR_f>64))=64; %to provent Nan numerical error
            LR_f(LR_f<-64)=-64; %to provent Nan numerical error

            LR_p = wb_LDPC_Decoder(H_Mesh,Q1,Q2,Num_iteration,LR_f);
         
            %Decision
            LDPCDecoderOut=(LR_p>0);            
            
            %compute BER only for message bits
            NumOfError=sum(sum(mod(LDPCDecoderOut(1:Ns)+MessageMatrix(index,:),2)))+NumOfError;
        end %end of index = 1:size(LR_fMatrix,1)
        
        %BER breaking criterion
        if (NumOfError>=1000)
            break;
        end %end of if
    end %end for repeat = 1:1e50
    
    BER = NumOfError/(repeat*noAveSamples*k*Ns);
    fprintf('BER = %f\n',BER);
    BERRecord = [BERRecord BER];
    save BERRecord BERRecord
    
    if (BER<=1e-4)
        break;
    end
end %end of EEsN02

figure();
semilogy(EbNodB_array,BERRecord);
grid on;
axis([0 8 1e-6 1e0]);
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');

% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix');
% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/LDPC_decoder');
%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');

% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');
% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns100_icr8_dv8_dc9');
% rmpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns200_icr12_dv6_dc13');


