%DAS(Distributed Assistant Antenna System) Network performance simulation
clc;clear;
%%
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
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
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns100_icr8_dv8_dc9');
% case 3
Ns = 200; Nr = 100; dv = 6; dc = 13;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns200_icr12_dv6_dc13');
code_rate = Ns / (Ns+Nr);
%%
%%%%# of inner iterations in LDPC decoder
Num_iteration=20; 

AttackPercent=0.05; %0.05 %# % of relay nodes are attached
Pinsist=0.5;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2

EsNodB_array = -4:1:6;
% EsNodB_array = 2:3;
% EbNodB_array = -2:1:5;

NoAttackedRelayNode = ceil(Nr*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
% tmp_candi = randperm(Nr);
% tmp_pos = tmp_candi(1:NoAttackedRelayNode);
% tmp_pos = sort(tmp_pos);
tmp_pos = 1:NoAttackedRelayNode;

BERRecord=[];
for idx=1:length(EsNodB_array)
%     EbNodB = EbNodB_array(idx); %unit in dB
%     EbNo = 10.^(EbNodB./10);
%     EsNo = EbNo * code_rate;
%     EsNodB = 10*log10(EsNo);    
%     No = 1/EsNo;
    
    EsNodB = EsNodB_array(idx); %unit in dB
    EsNo = 10.^(EsNodB./10);
    No = 1/EsNo;
        
    fprintf('Eb/No = %d [dB]\n',EsNodB);
    fprintf('Noise power = %f \n',No);

    NumOfError = 0;        
    
    for repeat = 1:1e6
        Message = randsrc(1, Ns, [0 1]);
%         Message = zeros(1, Ns);
        Codeword_origin = mod(Message*G_NC, 2);
        
        AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]);
        Codeword = Codeword_origin;
        Codeword(Ns + tmp_pos) = mod(Codeword_origin(Ns + tmp_pos) + AttackPatch, 2);
        
        noise = randn(1,Ns+Nr)*sqrt(No/2);
        signal = 2*Codeword - 1;
        y = signal + noise;
                
        LR_f = (4*EsNo).*y;
        %decoding process
        LR_f((LR_f>64))=64; %to provent Nan numerical error
        LR_f(LR_f<-64)=-64; %to provent Nan numerical error

        LR_p = wb_LDPC_Decoder(H_Mesh,Q1,Q2,Num_iteration,LR_f);
         
        %Decision
        LDPCDecoderOut = (LR_p>0);            
        NumOfError = sum(sum(mod(LDPCDecoderOut(1:Ns)+Message,2))) + NumOfError;
        
        %BER breaking criterion
        if (NumOfError>=1000)
            break;
        end %end of if
    end %end for repeat = 1:1e50
    
    BER = NumOfError/(repeat*Ns);
    fprintf('BER = %f\n',BER);
    BERRecord = [BERRecord BER];
    save BERRecord BERRecord
    
    if (BER<=1e-5)
        break;
    end
end %end of EEsN02

figure();
semilogy(EsNodB_array,BERRecord);
grid on;
% axis([0 8 1e-6 1e0]);
%%
% rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
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