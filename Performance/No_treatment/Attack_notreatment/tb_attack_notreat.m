%DAS(Distributed Assistant Antenna System) Network performance simulation
clc;clear;close all;
%%
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');
%%
% Ns = 50; Nr = 100; dv = 10; dc = 6;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
% addpath('/Users/woongbi/Dropbox/matlab_network_security/wblee_2018/H_matrix/Fixed_Nr100/Ns50_icr5_dv10_dc6');
Ns = 100; Nr = 100; dv = 8; dc = 9;
addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
% Ns = 200; Nr = 100; dv = 6; dc = 13;
% addpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');
code_rate = Ns / (Ns+Nr);

%System Parameters
% Ns=100; %# of Senders
% Nd=100; %# of distributed relays
% ICR=5;%# of incoming connection per relay 
% %Enter 1 for no time-domain coding
% n=1;  %# of output coded bits in time-domain
k=1;  %# of input bits in time-domain
% H_j=1;%# of 1's in a column in the H matrix
% H_k=1;%# of 1's in a row i the H matrix
%%%%# of inner iterations in LDPC decoder
Num_iteration = 20;

AttackPercent=0.05; %0.05 %# % of relay nodes are attached
Pinsist=1;        %Insistance of the attacker to flip the binary message out of all attacked relay nodes

load GG GG
load G_NC G_NC
load H_Mesh H_Mesh
load Q1 Q1
load Q2 Q2
%%%
EbNodB_array = 0:1:2;
% EbNodB_array = 30;
BERRecord=[];

NoAttackedRelayNode = ceil(Nr*AttackPercent); % No. of Attacked Relay Nodes among total No. of Relay Nodes
tmp_candi = randperm(Nr);
tmp_pos = tmp_candi(1:NoAttackedRelayNode);
tmp_pos = sort(tmp_pos);
tmp_pos = 1:NoAttackedRelayNode;

for idx=1:length(EbNodB_array)
    EbNodB = EbNodB_array(idx); %unit in dB
    EbNo = 10.^(EbNodB./10);
    EsNo = EbNo * code_rate;
    EsNodB = 10*log10(EsNo);    
    Es = 1/code_rate
    No = 1/EsNo;
    
    fprintf('Eb/No = %d [dB]\n',EbNodB);
    fprintf('Noise power = %f \n',No);
    NumOfError = 0;    
    
    for repeat = 1:1e6
        MessageMatrix = randsrc(1, Ns, [0 1]);
        Codeword = MessageMatrix;
        MeshCodeword=[];
        MeshCodewordTemp = mod(Codeword*G_NC ,2);
        AttackPatch = randsrc(1, NoAttackedRelayNode, [0 1; (1-Pinsist) Pinsist]); % Choose first No. of Attacked Relay Nodes
                
                %AttackedMeshCodeword
%                 MeshCodewordTemp(Ns+1:Ns+NoAttackedRelayNode) = mod(MeshCodewordTemp(Ns+1: Ns+NoAttackedRelayNode)+AttackPatch, 2);
        MeshCodeword = MeshCodewordTemp;
        MeshCodeword(Ns + tmp_pos) = mod(MeshCodeword(Ns + tmp_pos) + AttackPatch, 2);
        
        %Channel effect
%             EsN02=EEsNO2(IEsNO2); %unit in dB
%             yMatrix = awgn((2*MeshCodeword - 1), EsNO2);
        noise = sqrt(No/2)*randn(1,Ns+Nr);
        sig = (2*MeshCodeword - 1);
        y = sig + noise;
%         y = sig;
            
            %LDPC decodor
            %initial:        
            
%             EsNO2_lin = 10^(EsNO2/10);
%             LR_f = (2*EsNO2_lin).*y;
        LR_f = (4*EsNo)*y;
            
            %decoding process
        LR_f((LR_f>64))=64; %to provent Nan numerical error
        LR_f(LR_f<-64)=-64; %to provent Nan numerical error

        LR_p=wb_LDPC_Decoder(H_Mesh,Q1,Q2,Num_iteration,LR_f);
         
            %Decision
        LDPCDecoderOut=(LR_p>0);            
            %compute BER only for message bits            
        NumOfError=sum(sum(mod(LDPCDecoderOut(1:Ns)+MessageMatrix,2))) + NumOfError;
                
        %BER breaking criterion
        if (NumOfError>=1000)
            break;
        end %end of if
    end %end for repeat = 1:1e50
    
    BER = NumOfError/(repeat*Ns);
    fprintf('BER = %f\n',BER);

    BERRecord = [BERRecord BER];    
    save BERRecord BERRecord
    
    if (BER<=1e-4)
        break;
    end
end %end of EEsN02
%%
figure();
semilogy(EbNodB_array,BERRecord);
grid on;
axis([0 9 1e-6 1e0]);
%%
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\LDPC_decoder');

rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns50_icr5_dv10_dc6');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns100_icr8_dv8_dc9');
rmpath('C:\Users\wbleelab\Dropbox\matlab_network_security\wblee_2018\H_matrix\Fixed_Nr100\Ns200_icr12_dv6_dc13');